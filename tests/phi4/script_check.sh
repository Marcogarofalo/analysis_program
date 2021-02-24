#!/bin/bash

bad=()
cd jackknife
for f in `ls | grep jack | grep -v reference`
do 
   diff -q $f "$f"_reference > /dev/null || bad+=( "$f" )
done 
cd ..


cd out
for f in `ls | grep -v reference`
do 
   diff -q $f "$f"_reference_jack > /dev/null || bad+=( "$f" )
done
cd ..

if [[ -n "${bad-}" ]]; then
  echo -e "\nDifferences found in:\n"
  for path in "${bad[@]}"; do
    echo "  - $path"
  done

  exit 1
fi
