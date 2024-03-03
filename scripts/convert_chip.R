suppressPackageStartupMessages(library(tidyverse))
library( RSQLite )
library( argparse )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_args = function() {
	parser = ArgumentParser( description = 'Combine multiple tabulate-mismatches outputs into one file.' )
	parser$add_argument(
		"--manifest",
		type = "character",
		help = "Name of input .sqlite file to load manifest annotation from",
		default = "data/GAMCC/microarray/eurofins_pmra/manifest/TFS-Assets_LSG_Support-Files_Axiom_PMRA/Axiom_PMRA.na35.annot.db"
	)
	parser$add_argument(
		"--genotypes",
		type = "character",
		help = "Path to Eurofins genotype file",
		default = "data/GAMCC/microarray/eurofins_pmra/S1025_2023_01-02_Genotype.txt.gz"
	)
	parser$add_argument(
		"--samples",
		type = "character",
		help = "Path to Eurofins samples file",
		default = "data/GAMCC/microarray/eurofins_pmra/S1025_2023_01-02_Sample_Table.txt"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Name of output .vcf file",
		required = TRUE
	)
	return( parser$parse_args() )
}

args = parse_args()

# echo( "++ Loading manifest annotation from %s...\n", args$manifest )
annot.db = dbConnect( dbDriver( "SQLite" ), args$manifest ) 
annotation = as_tibble( dbGetQuery( annot.db, "SELECT * FROM Annotations" ))
# It turns out all PMRA SNPs are annotated on the + strand, which is helpful.
# However let's check this here:
stopifnot( length( which( annotation$Strand != '+' )) == 0)
# echo( "++ Ok, loaded %d annotation records.\n", nrow(annotation))
# echo(
# 	"++ %d records had Stop != Start, max length was %d...\n",
# 	length( which( annotation$Start != annotation$Stop )),
# 	max( annotation$Stop - annotation$Start, na.rm = T )
# )

# echo( "++ Loading genotypes from %s...\n", args$genotypes )
X = read_tsv( args$genotypes, comment = '#' )

metadata = X[,c(1,188:212)]
G = as.matrix( X[,2:187])
colnames(G) = gsub( ".CEL_call_code", "", colnames(G))

# echo( "++ Loading samples from %s...\n", args$samples )
samples = read_tsv( args$samples, comment = '#' )
samples$sample_id = gsub( ".CEL", "", samples[[1]], fixed = T )
M = match( colnames(G), samples$sample_id )
stopifnot( length( which( is.na( M ))) == 0 )
stopifnot( length( M ) == ncol(G) )
samples = samples[M,]

# echo( "++ Ok, loaded genotypes for %d samples and %d variants.\n", nrow(samples), nrow(G))

# Match up genotypes to annotation
# echo( "++ Matching genotypes to annotations...\n" )
M = match( metadata$probeset_id, annotation$ProbeSet_ID )
stopifnot( length( which( is.na(M))) == 0 )
# echo( "++ Ok, all variants matched.\n" )

matched.annotation = annotation[M,]

# We IGNORE a few records that are apparently for multi-allelic
# variants, but which have two records with the second expressed as two alt alleles.
# => hard to process sensibly without further work.
info = tibble( ProbeSet_ID = matched.annotation$ProbeSet_ID, result = NA )
info$result[wNoRef] = "no_ref_allele"

wNoRef = which( matched.annotation$Ref_Allele == '.' | is.na( matched.annotation$Ref_Allele ) )
# echo( "++ I will ignore %d variants with . ref allele, (usually means multiple alleles).\n", length( wNoRef ))


# Make sure all remaining variants are either ref/alt or alt/ref
{
	w = which(
		# Not ref/alt
		( matched.annotation$Ref_Allele == matched.annotation$Allele_A & matched.annotation$Alt_Allele == matched.annotation$Allele_B )
		# and not alt/ref
		| ( matched.annotation$Ref_Allele == matched.annotation$Allele_B & matched.annotation$Alt_Allele == matched.annotation$Allele_A )
	)
	stopifnot( length(w) == nrow(matched.annotation) - length(wNoRef))
}

allele.sets = unique( matched.annotation[, c( "Ref_Allele", "Alt_Allele")] )
# print( allele.sets )

# echo( "++ Translating genotype calls to VCF format...\n" )
translated = G
translated[,] = "./."
counts = c( total = 0, fail = 0, noref = length( wNoRef ) )

for( i in 1:nrow( allele.sets )) {
	set = allele.sets[i,]
	# echo( "  --: %s...\n", paste( set, collapse = " > " ))

	translation = c()
	translation[sprintf( "%s/%s", set$Ref_Allele, set$Ref_Allele )] = '0/0'
	translation[sprintf( "%s/%s", set$Ref_Allele, set$Alt_Allele )] = '0/1'
	translation[sprintf( "%s/%s", set$Alt_Allele, set$Ref_Allele )] = '0/1'
	translation[sprintf( "%s/%s", set$Alt_Allele, set$Alt_Allele )] = '1/1'
	translation[ "---"] = "./."

	w = which(
		!is.na( matched.annotation$Ref_Allele )
		& matched.annotation$Ref_Allele == set$Ref_Allele
		& matched.annotation$Alt_Allele == set$Alt_Allele
	)
	for( j in w ) {
		translated[j,] = sapply(
			G[j,],
			function(g) {
				translation[g]
			}
		)
	}
	info$result[w] = "ok"
	{
		fails = unique( which( is.na( translated[w,,drop=F] ), arr.in = TRUE )[,1])
		info$result[fails] = "fail"

		if( length(fails) > 0 ) {
			echo( "!! there were %d failures of %d variants (%.1f%%).\n", length(fails), length(w), 100 * length(fails) / length(w))
		}
		translated[w,][is.na(translated[w,])] = "./."
		counts['total'] = counts['total'] + length(w)
		counts['fail'] = counts['fail'] + length(fails)
	}
	echo(
		"++ ...processed %d variants so far, of which %d (%.1f%%) had annotation-data mismatches.\n",
		counts['total'],
		counts['fail'],
		100 * counts['fail'] / counts['total']
	)
}

# echo(
# 	"++ Ok, counts were: %d total, %d fail, %d no reference allele.\n",
# 	counts['total'],
# 	counts['fail'],
# 	counts['noref']
# )

# echo( "++ Ok, writing info output...\n" )
info$ref_allele = matched.annotation$Ref_Allele
info$alt_allele = matched.annotation$Alt_Allele
write_tsv( info, file = gsub( '[.]vcf', '.info.tsv', args$output, ))

# echo( "++ Ok, forming VCF output...\n" )
# We remove missing-ref-allele variants here
# As I don't think there can be a missing ref allele in VCF output.
if( length(wNoRef) > 0 ) {
	G = G[-wNoRef,]
	matched.annotation = matched.annotation[-wNoRef,]
	metadata = metadata[-wNoRef,]
}

vcf.data = tibble(
	CHROM = matched.annotation$Chr_id,
	POS = matched.annotation$Start,
	ID = matched.annotation$ProbeSet_ID,
	REF = matched.annotation$Ref_Allele,
	ALT = matched.annotation$Alt_Allele,
	QUAL = ".",
	FILTER = ".",
	INFO = sprintf(
		"PID=%s;PAR=%s;FLANK=%s",
		matched.annotation$ProbeSet_ID,
		ifelse(
			matched.annotation$Chr_id == 24,
			sprintf( "%d", matched.annotation$ChrX_PAR ),
			"."
		),
		matched.annotation$Flank
	),
	FORMAT = "GT"
)
vcf.data$ID[ !is.na( matched.annotation$dbSNP_RS_ID )] = matched.annotation$dbSNP_RS_ID[ !is.na( matched.annotation$dbSNP_RS_ID )]
table( vcf.data$CHROM, useNA="always")

vcf.data$CHROM = factor(
	vcf.data$CHROM,
	levels = 1:26
)
levels( vcf.data$CHROM ) = c( 1:22, NA, "X", "Y", "MT" )

ordering = order( vcf.data$CHROM, vcf.data$POS, vcf.data$REF, vcf.data$ALT )

# echo( "++ Ok, writing VCF output to %s...\n", args$output )
outputConnection = file( args$output )
writeLines(
	c(
		"##fileformat=VCFv4.2",
		'##INFO=<ID=PID,Type=String,Number=1,Description="Probeset ID from Affymetrix annotation file">',
		'##INFO=<ID=PAR,Type=String,Number=1,Description="For chromosome X variants, is the variant annotated in the PAR region?">',
		'##INFO=<ID=FLANK,Type=String,Number=1,Description="Flanking sequence as given in the annotatino file.">',
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
		sprintf(
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT%s",
			paste( sprintf( "\t%s", colnames(translated) ), collapse = "" )
		)
	),
	con = outputConnection
)

write_tsv(
	bind_cols(
		vcf.data[ordering,], translated[ordering,]
	),
	file = outputConnection,
	col_names = FALSE,
	append = TRUE
)

samples_output = gsub( ".vcf", ".sample", args$output, fixed = TRUE )
# echo( "++ Ok, writing samples output to %s...\n", samples_output )
writeLines(
	c(
		"ID\tstatus\tDQC\tQC_call_rate\tcall_rate\tQC_het_rate\thet_rate\tQC_computed_gender\tplate\twell",
		"0\tD\tC\tC\tC\tC\tC\tD\tD\tD"
	),
	con = samples_output
)
write_tsv(
	samples[ , c(
		"sample_id",
		"Pass/Fail", "DQC", "QC call_rate", "call_rate", "QC het_rate", "het_rate",
		"QC computed_gender",
		"affymetrix-plate-barcode", "affymetrix-plate-peg-wellposition"
	)],
	file = samples_output,
	col_names = FALSE,
	append = TRUE
)
# echo( "++ Success.\n" )
# echo( "++ Thanks for using convert.R\n" )

