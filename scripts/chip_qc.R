library( dbplyr )
library( dplyr )
library( RSQLite )
db = dbConnect( dbDriver( "SQLite" ), "results/microarray/GRCh38/qc/S1025_2023_01-02_Genotype.qc.sqlite" )
snp_stats = ( db %>% tbl( "autosomes" ) %>% collect() )
sample_stats = ( db %>% tbl( "sample_stats" ) %>% collect() )
pcs = ( db %>% tbl( "PCsView" ) %>% filter( analysis %in% c( 'PCs:thin_1bp:exclude-duplicates', 'PCs:thin_1bp:all' )) %>% collect() )
kinship = read_tsv( "results/microarray/GRCh38/qc/PCs/S1025_2023_01-02_kinship_thin_1bp.exclude-duplicates.tsv.gz", comment = '#' )

pdf( file = "tmp/images/PDRA_GRCh38_stats.pdf" )
hist( snp_stats$missing_proportion, 20, xlab = "Variants: propn missing" )
hist( snp_stats$minor_allele_frequency, 50, xlab = "Variants: MAF" )
hist( -log10( snp_stats$HW_exact_p_value), 50, xlab = "Variants: -log10 HW P-value" )
hist( sample_stats$missing_proportion, 20, xlab = "Samples: propn missing genotypes" )
hist( sample_stats$heterozygous_proportion, 20, xlab = "Samples: propn het genotypes" )
( 
	ggplot( data = sample_stats )
	+ geom_point( aes( x = missing_proportion, y = heterozygous_proportion ))
)
( 
	ggplot( data = pcs %>% filter( analysis == 'PCs:thin_1bp:all' ))
	+ geom_point( aes( x = PC_1, y = PC_2, colour = ))
	+ ggtitle( "PCs across all samples")
	+ theme_minimal()
)
( 
	ggplot( data = pcs %>% filter( analysis == 'PCs:thin_1bp:exclude-duplicates' ))
	+ geom_point( aes( x = PC_1, y = PC_2, colour = ))
	+ ggtitle( "PCs after removing two apparent duplicates")
	+ theme_minimal()
)
pairs(
	( pcs %>% filter( analysis == 'PCs:thin_1bp:exclude-duplicates' ))[,sprintf( "PC_%d", 1:5 )],
	pch = 19
)
hist( (kinship %>% filter( sample_1 == sample_2 ))$value, 20, xlab = "Estimated kinship", main = "Kinship with self" )
hist( (kinship %>% filter( sample_1 != sample_2 ))$value, 50, xlab = "Estimated kinship", main = "Kinship with other" )

dev.off()
