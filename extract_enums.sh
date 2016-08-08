#!/bin/bash
inone=0
words=""
enums=""
echo "#include <enum_str.h>" > enum_str.C
while read line ; do
    
    #echo "$line"
    if [[ $inone -eq 1 ]] ; then

        regex='}'
        if [[ $line =~ $regex ]] ; then
            inone=0
            echo "const char* ${enum}_str[] = {"
            for w in $words ; do
                echo '"'$w'"',
            done
            echo '""};'
            enums="${enums} $enum"
            words=""
            continue
        fi
        words="${words} "`echo "$line" | sed -e "s:,: :g" | sed -e "s:/.*::g"`

    fi

    regex='^enum[[:space:]]+[^[:space:]]+[[:space:]]*[{]'
    if [[ $line =~ $regex ]] ; then
        inone=1
        set -- $line
        enum=$2
    fi
done >> enum_str.C

echo "#pragma once" > enum_str.h
for w in $enums ;  do
    echo "extern const char* ${w}_str[];" >> enum_str.h
done


