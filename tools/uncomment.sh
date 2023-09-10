#!/bin/bash

    input_file="${MS}/valgrind-python-to-edit.supp"
    output_file="${MS}/valgrind-python.supp"
    inside_block=false

    rm -f "$output_file"

    uncomment_block=""
    old_block=""
    while IFS= read -r line; do
      stripped_line=$(echo "$line" | xargs -0)  # Remove leading/trailing white space

      uncomment_block="$uncomment_block\n${line:3}"
      old_block="$old_block\n$line"

      if [[ "$stripped_line" == "{" ]]; then
        inside_block=false
      fi

      if [[ "$stripped_line" == *"_PyObject_Free"* || "$stripped_line" == *"_PyObject_Realloc"* ]]; then
        inside_block=true
      fi

      if [[ "$stripped_line" == *"}" ]]; then
         if [ "$inside_block" == true ]; then
            echo -e "$uncomment_block" >> "$output_file"
         else
            echo -e "$old_block" >> "$output_file"
         fi
         uncomment_block=""
         old_block=""

      fi


    rm -f "$input_file"

    done < "$input_file"