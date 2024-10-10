
if [ $# -ne 3 ]; then
    echo "Error: Please enter 2input file name.and at least 1 outputfilename "
    exit 1
fi

cat "$1" > "$3"

cat "$2" >> "$3"

echo "Merged File is: $3"
cat "$3"

exit
