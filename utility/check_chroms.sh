source ${WORK_BGFS}/scripts/master/refs.env

if echo "$1" | grep -q "\."
then
    prefix=$(echo $1 | awk 'BEGIN{FS=OFS="."} NF--')
else
    prefix="$1"
fi

out="${prefix}.check_chroms"

printf "" > $out
while read l
do
    echo "$l = `grep -c -P "^$l" $1`" >> $out
done < $CRHOMS