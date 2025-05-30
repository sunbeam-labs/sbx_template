#!/usr/bin/env bash
while getopts a:n:u:d: flag
do
    case "${flag}" in
        a) author=${OPTARG};;
        n) name=${OPTARG};;
        u) urlname=${OPTARG};;
        d) description=${OPTARG};;
    esac
done

full_name="$author/$name"
echo "Author: $author";
echo "Project Name: $name";
echo "Project Full Name: $full_name";
echo "Project URL name: $urlname";
echo "Description: $description";

echo "Renaming project..."

original_full_name="sunbeam-labs/sbx_template"
original_name="sbx_template"
#original_urlname="project_urlname"
#original_description="project_description"
# for filename in $(find . -name "*.*") 
for filename in $(git ls-files) 
do
    if [[ $filename == .github/workflows/*.yml ]]; then
        continue
    fi
    sed -i "s/$original_full_name/$full_name/g" $filename
    sed -i "s/$original_name/$name/g" $filename
    #sed -i "s/$original_urlname/$urlname/g" $filename
    #sed -i "s/$original_description/$description/g" $filename
    echo "Renamed $filename"
done

mv sbx_template.smk $name.smk
mv envs/sbx_template_env.yml "envs/${name}_env.yml"

# This command runs only once on GHA!
rm -f .github/template.yml
rm -f .github/rename_project.sh
rm -f .github/workflows/rename_extension.yml
echo "Project renamed successfully!"