#!/usr/bin/awk -f
BEGIN {
  prev = ""
  prevline = ""
}

{
    if ($0 ~ "^@") {
    print $0
    next
    }
    if ($1 == prev) {
    print prevline
    print $0
    } else {
    prev = $1
    prevline = $0
    }
}