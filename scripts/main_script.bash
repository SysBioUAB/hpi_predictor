#!/bin/bash

                            # AQUÍ UN CUADRO INFORMATIVO CON EL NOMBRE DEL AUTOR #


                            # CAMBIAR PRESENTACIÓN Y NOMBRES DE VARIABLES #

	# Variables #

	## directories ##
	script_dir=$(pwd)
    base_dir=$(sed 's/scripts//' <<< $script_dir)
   	result_dir="${base_dir}"results/
   	host_dir="${base_dir}"datasets/hosts/
   	pathogen_dir="${base_dir}"datasets/pathogens/
   	pos_inter_dir="${base_dir}"datasets/positive_interactions/
	helper_dir_list=("$host_dir" "$pathogen_dir" "$pos_inter_dir")
   	host_list="${base_dir}"datasets/hosts/host_list
   	pathogen_list="${base_dir}"datasets/pathogens/pathogen_list
   	pos_dataset="${base_dir}"datasets/positive_interactions/phi_data.tsv 
   	def_descriptors_list="${base_dir}"datasets/descriptors/default_indices
   	custom_descriptors_list="${base_dir}"datasets/descriptors/custom_indices
   	av_tax="${base_dir}"datasets/available_taxon_IDs.tsv
   	helper_list=("$host_list" "$pathogen_list" "$av_tax")
   	desc_list="${base_dir}"datasets/descriptors/descriptors_list.tsv
	helper_names=("host" "pathogen")
	databases_dir="${base_dir}"databases/
	dualseqdb_file="${databases_dir}"dualseqdb_v1.zip
	dualseqdb_file_unz="${databases_dir}"dualseqdb_v1.tsv
	bacfitbase_file="${databases_dir}"bacfitbase_v1.zip
	bacfitbase_file_unz="${databases_dir}"bacfitbase_v1.tsv
	interactome_file="${result_dir}"temp_file
	phi_base_consensus_file="${databases_dir}"phi-base_consensus_phenotypes.fasta.psq
    interactome_subdirs="${result_dir}"temp_file_subdir
	declare -a organisms_user=()

# Ask user to choose host and pathogen # 

for i in {0..1}; do

    organisms=$(cat ${helper_list[i]})  
    organism_name=$(echo -e "$organisms\nOther" | zenity --title="List of available "${helper_names[i]}"  " --text="Please choose one from the list below. Choose other for custom organism" --list --column="${helper_names[i]}" 2>/dev/null)


    if [ -z "$organism_name" ] ; then exit; fi



# Enters here if an organism is chosen

    if [ "$organism_name" != "Other" ]; then
            organisms_user+=("$organism_name")
            if [ "$i" == 0 ] && [[ ! "$organism_name" =~ .*custom$.* ]] ; then
                host_taxonID=$(sed 's/^.\{1,\} \([[:digit:]]\{1,\}\)$/\1/' <<< "$organism_name")
#                echo "host taxon ID $host_taxonID"
            fi
# If Other is chosen as organism, give chance to choose between entering taxon ID or give custom file
   

    elif [ "$organism_name" == "Other" ]; then

        choice=$(zenity --title="Choosing custom organism" --text="Choose -Taxon ID- to download proteome from available proteomes in Uniprot. Choose -Custom file- to enter custom proteome file (File format example given)" --list "Taxon ID" "Custom file"  --column="Custom method" 2>/dev/null)

        if [ -z "$choice" ] ; then exit; fi

        if [ "$choice" == "Taxon ID" ]; then

            taxonID=$(zenity --entry --title="Enter taxon ID for ${helper_names[i]}" 2>/dev/null)
            if [ -z "$taxonID" ] ; then exit; fi
            var=$(grep -w "^$taxonID" ${helper_list[2]} | head -n 1)
            while [ -z "$var" ]; do
                taxonID=$(zenity --entry --title="Enter taxon ID for ${helper_names[i]}" --text="Last Taxon Id not available. Please check available Taxon ID list and try again:" 2>/dev/null)
                if [ -z "$taxonID" ] ; then exit; fi
                var=$(grep -w "^$taxonID" ${helper_list[2]})
            done
            organism_name=$(grep -w "^$taxonID" ${helper_list[2]} | head -n 1 | awk -F'\t' '{print $2}')

            if [ "$i" == 0 ]; then

                host_taxonID=$(echo "$taxonID")            
            fi

            organism_name2=$(echo "$organism_name" | sed  's/(.\+$//')
            organism_name3=$(echo "$organism_name2$taxonID")
 #           echo -e "$organism_name3"
            organisms_user+=("$organism_name3")
            var=$(grep "^$organism_name3$" ${helper_list[i]})

            if [ -z "$var" ]; then

                echo -e "Downloading proteome for $organism_name3 . This may take a while"
                organism_help=$(sed 's/ /_/g;s/)$//g;s/_(/_/;s/_\/_/_/g;s/\.$//;s/).(/_/g;s/\.//g;s/://g;s/\//_/g' <<< $organism_name3)
                organism_file=$(sed 's/ /_/g;s/$/_proteome.tsv/' <<< $organism_help)
                new_organism_dir=$(echo "${helper_dir_list[i]}$organism_help/")
               
                if [ ! -d "$new_organism_dir" ]; then mkdir "$new_organism_dir"; fi

                new_organism_file_fullpath=$(echo "$new_organism_dir$organism_file")
                curl -s "https://www.uniprot.org/uniprot/?query=organism%3A%22"$taxonID"%22&sort=score&columns=id,protein%20,organism,reviewed,taxon,length,sequence&format=tab" > "$new_organism_file_fullpath"
                reviewed=$(grep -w "reviewed" "$new_organism_file_fullpath" | wc -l)
                unreviewed=$(grep -w "unreviewed" "$new_organism_file_fullpath" | wc -l)
                all=$(bc <<< "scale=2; $reviewed+$unreviewed")
                echo -e "$reviewed reviewed proteins for $organism_name3"
                echo -e "$unreviewed unreviewed proteins for  $organism_name3"
            
                if [ "$reviewed" -eq 0 ]; then
   
                    zenity --warning --title="Warning" --text="<b> - $organism_name3 - </b> has no reviewed entries available. Unreviewed entries will be used instead" 2>/dev/null

                    awk -F'\t' -v OFS='\t' '{print $1, $2, $4, $5, $6}' "$new_organism_file_fullpath" > helper && rm "$new_organism_file_fullpath" && mv helper "$new_organism_file_fullpath"
                    new_organism_desc_dir=$(echo $new_organism_dir"descriptors")

                    if [ ! -d "$new_organism_desc_dir" ]; then mkdir "$new_organism_desc_dir"; fi
                else
                    prot=$(zenity --title="Proteome status selection" --text="Choose REVIEWED entries OR ALL entries (reviewed + unreviewed) for <b> $organism_name3 </b>" --list --separator " " --column="Status" --column="Number of entries" REVIEWED  $reviewed ALL $all  2>/dev/null)

                    if [ -z "$prot" ] ; then exit; fi

                    if  [ "$prot" == "ALL" ] ; then
                        awk -F'\t' -v OFS='\t' '{print $1, $2, $4, $5, $6}' "$new_organism_file_fullpath" > helper && rm "$new_organism_file_fullpath" && mv helper "$new_organism_file_fullpath"
                        new_organism_desc_dir=$(echo $new_organism_dir"descriptors")
                        if [ ! -d "$new_organism_desc_dir" ]; then mkdir "$new_organism_desc_dir"; fi
                    elif [ "$prot" == "REVIEWED" ] ; then
                        sed -e '1p' -e '/\treviewed/!d' "$new_organism_file_fullpath" | awk -F'\t' -v OFS='\t' '{print $1, $2, $4, $5, $6}' > helper && rm "$new_organism_file_fullpath" && mv helper "$new_organism_file_fullpath"
                        new_organism_desc_dir=$(echo $new_organism_dir"descriptors")
                        if [ ! -d "$new_organism_desc_dir" ]; then mkdir "$new_organism_desc_dir"; fi
                    fi
                fi


                if  ! grep -Fxq "$organism_name3" "${helper_list[i]}" ; then
                    echo "$organism_name3" >> "${helper_list[i]}"
                fi

            else

                echo -e  "$organism_name already in our database"

            fi 

        elif [ "$choice" == "Custom file" ]; then

            organism_name=$(zenity --entry --title="Custom host name" --text="Please enter host name separated by spaces. Last name must be custom" 2>/dev/null)
            if [ -z "$organism_name" ] ; then exit; fi


            if [[ ! "$organism_name" =~ .*custom$.* ]]; then exit; fi

            path=$(zenity --file-selection title="Custom proteome" --text="Please choose a valid file for your custom proteome" 2>/dev/null)
            if [ -z "$path" ] ; then exit; fi

            var=$(grep $'Entry\tOrganism\tOrganism ID\tLength\tSequence' $path 2>/dev/null)

            while [ "$var" != $'Entry\tOrganism\tOrganism ID\tLength\tSequence' ] ; do

                    path=$(zenity --file-selection title="Custom proteome" --text="Please choose a valid file for your custom proteome" 2>/dev/null)
                    if [ -z "$path" ] ; then exit; fi
                    var=$(grep $'Entry\tOrganism\tOrganism ID\tLength\tSequence' $path 2>/dev/null)

            done

            organism_help=$(sed 's/ /_/g;s/)$//g;s/_(/_/;s/_\/_/_/g;s/\.$//;s/).(/_/g;s/\.//g;s/://g' <<< $organism_name)
            organism_file=$(sed 's/ /_/g;s/$/_proteome.tsv/' <<< $organism_help)
            new_organism_dir=$(echo "${helper_dir_list[i]}$organism_help/")
            if [ ! -d "$new_organism_dir" ]; then mkdir "$new_organism_dir"; fi
            new_organism_file_fullpath=$(echo "$new_organism_dir$organism_file")
            cp "$path" "$new_organism_file_fullpath"
            new_organism_desc_dir=$(echo $new_organism_dir"descriptors")
            if [ ! -d "$new_organism_desc_dir" ]; then mkdir "$new_organism_desc_dir"; fi
            organisms_user+=("$organism_name")

            if  ! grep -Fxq "$organism_name" "${helper_list[i]}" ; then
                echo "$organism_name" >> "${helper_list[i]}"
            fi

        fi
   
    fi

done 

# echo "host ${organisms_user[0]}"
# echo "pathogen ${organisms_user[1]}"

# echo "${helper_dir_list[0]}" 
# echo "${helper_dir_list[1]}"

#### 3. Ask user to choose descriptors

answer=$(zenity --title="Physicochemical descriptors" --text="Please choose between default or custom descriptors" --list "Default" "Custom" --column="Descriptor mode" 2>/dev/null) 
if [ -z "$answer" ] ; then exit; fi


if [ "$answer" == "Default" ]; then

    desc=$(echo "default")

elif [ "$answer" == "Custom" ]; then

    descriptors=$(cat "$desc_list" | zenity --title="Choosing custom descriptors" --list --multiple --column="Descriptor" 2>/dev/null)
    if [ -z "$descriptors" ] ; then exit; fi
    echo "$descriptors" | sed 's/|/\n/g' | awk '{print $1}' > $custom_descriptors_list

fi  

lagmax=200
cutoff_p="0.4"

false_positive=$(zenity --list --title="False Positive Rate" --text="Please choose the percentage of false positives allowed to be part of predicted interactions" "0.1" "0.01" "0.001" "0.0001" "0.00001" --column="FPR" 2>/dev/null)
if [ -z "$false_positive" ] ; then exit; fi

consensus_predictors=$(zenity --list --title="Minimum number of of agreeing descriptors" --text="Please choose the percentage of Physicochemical descriptors that need to agree on a prediction in order to allow it to be part of the predicted interactome" "50" "60" "70" "80" "90" "100" --column="% of descriptors" 2>/dev/null)
if [ -z "$consensus_predictors" ] ; then exit; fi
# percentage of predictors agreeing on same prediction to be considered as a predicted interaction

consensus_predictors_p=$(bc <<< "scale=2; $consensus_predictors/100")

### invoking PROCESSING SCRIPT
#echo "${organisms_user[0]}"
#echo "${organisms_user[1]}"
#echo "${helper_dir_list[0]}"
#echo "${helper_dir_list[1]}"
#echo "${helper_dir_list[2]}"
#echo "$lagmax"
#echo "$cutoff_p"
#echo "$result_dir"
#echo "$false_positive"
#echo "$consensus_predictors_p"
#echo "$custom_descriptors_dir"
#echo "$taxonID"

time Rscript interaction_prediction_script.R "${organisms_user[0]}" "${organisms_user[1]}" "${helper_dir_list[0]}" "${helper_dir_list[1]}" "${helper_dir_list[2]}" "$lagmax" "$cutoff_p" "$result_dir" "$base_dir" "$pos_dataset" "$def_descriptors_list" "$custom_descriptors_list" "$false_positive" "$consensus_predictors_p"

interactome=$(cat $interactome_file)
if [ "$interactome" == "No consensus interactome generated" ]; then echo "No entries generated for consensus interactome. End of analysis." && rm "$interactome_file" && exit; fi

interactome_length=$(cat "$interactome" | wc -l)
if [ "$interactome_length" -le 6 ]; then echo "Not enough entries generated. End of the analysis." && rm "$interactome_file" && exit; fi
#if [ "$interactome_length" -eq 1 ]; then echo "No entries generated for consensus interactome" && rm "$interactome" && exit; fi
intermediate_files=$(cat $interactome_subdirs | sed -n 2p)
interactome_sub_files=$(cat $interactome_subdirs | sed -n 1p)

#uncompress databases

cd  $databases_dir


host=$(sed 's/ /_/g' <<< ${organisms_user[0]})
pathogen=$(sed 's/ /_/g' <<< ${organisms_user[1]})

if [ ! -f "$dualseqdb_file_unz" ] ; then
    unzip "$dualseqdb_file"
    #disccard entries with no sequences and isoforms (differente sequences, because of the isoforms, have same values associated, makes no sense)

    awk -F'\t' -v OFS='\t' '{ if ($19 != "[NULL]") { print } }'  dualseqdb_v1.tsv | awk -F'\t' -v OFS='\t'  '!seen[$1, $14]++' > dualseqdb_v1_withseq.tsv
    awk '{ if ($2 == "pathogen") { print } }' dualseqdb_v1_withseq.tsv | awk -F'\t' -v OFS='\t' '{print $9, $16, $17, $19}' | awk '!seen[$4]++' > filtered_dualseqdb_v1_withseq_bac_with_scores.tsv
    awk -F'\t' -v OFS='\t' '{print $1, $4}' filtered_dualseqdb_v1_withseq_bac_with_scores.tsv | sed 's/^/>/' | uniq  | tr '\t' '\n' > filtered_dualseqdb_v1_withseq_bac.fasta
    awk '{ if ($2 == "pathogen") { print } }' dualseqdb_v1_withseq.tsv | awk -F'\t' -v OFS='\t' '{print $9, $16, $17, $19}' > filtered_dualseqdb_v1_withseq_bac_with_scores.tsv

    awk '{ if ($2 != "pathogen") { print } }' dualseqdb_v1_withseq.tsv | awk -F'\t' -v OFS='\t' '{print $9, $16, $17, $19}' | awk '!seen[$4]++' | sed -n '2,$'p > filtered_dualseqdb_v1_withseq_human_with_scores.tsv
    awk -F'\t' -v OFS='\t' '{print $1, $4}' filtered_dualseqdb_v1_withseq_human_with_scores.tsv | sed 's/^/>/' | tr '\t' '\n' > filtered_dualseqdb_v1_withseq_human.fasta
    awk '{ if ($2 != "pathogen") { print } }' dualseqdb_v1_withseq.tsv | awk -F'\t' -v OFS='\t' '{print $9, $16, $17, $19}' | sed -n '2,$'p > filtered_dualseqdb_v1_withseq_human_with_scores.tsv
    makeblastdb -in filtered_dualseqdb_v1_withseq_bac.fasta -dbtype prot
    makeblastdb -in filtered_dualseqdb_v1_withseq_human.fasta -dbtype prot
fi

    falsepos_var=$(sed 's/\./_/g' <<< $false_positive)

    pred_interactome_by_some_pred_bac=$(echo "$intermediate_files"bacteria.fasta)
    pred_interactome_by_some_pred_host=$(echo "$intermediate_files"host.fasta)

#    pred_interactome_by_some_pred_bac=$(echo "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_bac.tsv)
#    pred_interactome_by_some_pred_host=$(echo "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_host.tsv)

    if [ ! -f "$pred_interactome_by_some_pred_bac" ]; then
        awk -F'\t' -v OFS='\t'  '{print $2, $7}' "$interactome" | sed -n '2,$'p | sed 's/^/>/' | sort | uniq | sort | tr '\t' '\n' > "$pred_interactome_by_some_pred_bac"

    fi
    if [ ! -f "$pred_interactome_by_some_pred_host" ]; then

        awk -F'\t' -v OFS='\t'  '{print $1, $23}' "$interactome" | sed -n '2,$'p | sed 's/^/>/' | sort | uniq | sort | tr '\t' '\n' > "$pred_interactome_by_some_pred_host"

    fi


########## bacfitbase ##########

if [ ! -f "$bacfitbase_file_unz" ] ; then
    unzip "$bacfitbase_file"
    awk -F'\t' -v OFS='\t' '{print $9, $19 }' bacfitbase_v1.tsv | sed -n '2,$'p | sed 's/^/>/' | uniq | tr '\t' '\n' > bacfitbase.fasta
    #create blastp database out of bacfitbase database
    makeblastdb -in bacfitbase.fasta -dbtype prot

fi

#create blaspt database out of phi-base

if [ ! -f "$phi_base_consensus_file" ] ; then
    makeblastdb -in phi-base_consensus_phenotypes.fasta -dbtype prot
fi

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#blast PAO1 pathogen fasta against bacfitbase fasta:


pred_interactome_by_some_pred_bac_blast=$(echo "$intermediate_files"bacteria_vs_bacfitbase.blast)
#pred_interactome_by_some_pred_bac_blast=$(echo "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_bac_blast.tsv)

#echo "$pred_interactome_by_some_pred_bac"

    if [ ! -f "$pred_interactome_by_some_pred_bac_blast" ]; then
        blastp -query $pred_interactome_by_some_pred_bac  -db bacfitbase.fasta -evalue 1e-50  -outfmt 6 | awk -F'\t' -v OFS='\t'  '{print $1, $2, $3, $11, $12}' > "$pred_interactome_by_some_pred_bac_blast"

    
# combine combine interactome and bacfitbase database

        query=$(awk -F'\t' -v OFS='\t' '{print $1}' "$pred_interactome_by_some_pred_bac_blast")
        awk '{print $2}' "$pred_interactome_by_some_pred_bac_blast" > template
        grep -wf template bacfitbase_v1.tsv | awk -F'\t' -v OFS='\t' '{print $9, $16, $17}' > templates
        templates=$(cat templates)
        echo "$templates" | while read p; do
            templates_query=$(echo "$p" | cut -f1)
            templates_fitness=$(echo "$p" | cut -f2)
            templates_pval=$(echo "$p" | cut -f3)
            helper=$(grep "$templates_query" "$pred_interactome_by_some_pred_bac_blast")
            freq=$(echo "$helper" | wc -l)
            templates_fitness2=$(for i in $(seq 1 $freq); do echo "$templates_fitness"; done)
            templates_pval2=$(for i in $(seq 1 $freq); do echo "$templates_pval"; done)
            paste <(echo "$helper") <(echo "$templates_fitness2") <(echo "$templates_pval2");
        done >"$intermediate_files"bacteria_vs_bacfitbase_values.tsv
    fi



# blast PAO1 pathogen fasta against phi-base database

    if [ ! -f phi-base_current.csv ]; then

        curl -s "https://raw.githubusercontent.com/PHI-base/data/master/releases/phi-base_current.csv" > phi-base_current.csv
    fi

    if [ ! -f phi-base_current.fasta ]; then
    
         curl -s "https://raw.githubusercontent.com/PHI-base/data/master/releases/phi-base_current.fas" > phi-base_current.fasta
    fi


    if [ ! -f phi-base_consensus_phenotypes.tsv ]; then

        Rscript "$script_dir"filter_phi-base.R phi-base_current.csv phi-base_consensus_phenotypes.tsv
    fi

    if [ ! -f phi-base_consensus_phenotypes.fasta ]; then
   
        helper=$(awk -F'\t' -v OFS='\t' '{print $1}' phi-base_consensus_phenotypes.tsv | sed -n '2,$'p)
        echo "$host_list" | while read p; do grep -A1 "$p" phi-base_current.fasta; done > phi-base_consensus_phenotypes.fasta

    fi

    pred_interactome_by_some_pred_bac_phi_blast=$(echo "$intermediate_files"bacteria_vs_phi_base.blast)
#    pred_interactome_by_some_pred_bac_phi_blast=$(echo "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_phi_base.tsv)
    if [ ! -f "$pred_interactome_by_some_pred_bac_phi_blast" ]; then
        blastp -query $pred_interactome_by_some_pred_bac -db phi-base_consensus_phenotypes.fasta -evalue 1e-50  -outfmt 6 | awk -F'\t' -v OFS='\t'  '{print $1, $2, $3, $11, $12}' > "$pred_interactome_by_some_pred_bac_phi_blast"

        query=$(awk -F'\t' -v OFS='\t' '{print $1}' "$pred_interactome_by_some_pred_bac_phi_blast")
        awk '{print $2}' "$pred_interactome_by_some_pred_bac_phi_blast" > template
        grep -wf template phi-base_consensus_phenotypes.tsv | awk -F'\t' -v OFS='\t' '{print $1, $2, $3}' > templates
        templates=$(cat templates)
        echo "$templates" | while read p; do
            templates_query=$(echo "$p" | cut -f1)
            templates_tag=$(echo "$p" | cut -f2)
            templates_score=$(echo "$p" | cut -f3)
            helper=$(grep "$templates_query" "$pred_interactome_by_some_pred_bac_phi_blast")
            freq=$(echo "$helper" | wc -l)
            templates_tag2=$(for i in $(seq 1 $freq); do echo "$templates_tag"; done)
            templates_score2=$(for i in $(seq 1 $freq); do echo "$templates_score"; done)
            paste <(echo "$helper") <(echo "$templates_tag2") <(echo "$templates_score2");
        done > "$intermediate_files"bacteria_vs_phi_base_values.tsv
    fi


# blast PAO1 pathogen fasta and human fasta against dualseqDB database


pred_interactome_by_some_pred_bac_dualseq_blast=$(echo "$intermediate_files"bacteria_vs_dualseq.blast)
pred_interactome_by_some_pred_host_dualseq_blast=$(echo "$intermediate_files"host_vs_dualseq.blast)

#pred_interactome_by_some_pred_bac_dualseq_blast=$(echo "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_bac_dualseq_blast.tsv)
#pred_interactome_by_some_pred_host_dualseq_blast=$(echo "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_host_dualseq_blast.tsv)


if [ ! -f "$pred_interactome_by_some_pred_bac_dualseq_blast" ]; then
    blastp -query $pred_interactome_by_some_pred_bac -db filtered_dualseqdb_v1_withseq_bac.fasta -evalue 1e-50  -outfmt 6 | awk -F'\t' -v OFS='\t'  '{print $1, $2, $3, $11, $12}' > "$pred_interactome_by_some_pred_bac_dualseq_blast"

    #bacteria
    query=$(awk -F'\t' -v OFS='\t' '{print $1}' "$pred_interactome_by_some_pred_bac_dualseq_blast")
    awk '{print $2}' "$pred_interactome_by_some_pred_bac_dualseq_blast" > template
    grep -wf template filtered_dualseqdb_v1_withseq_bac_with_scores.tsv | awk -F'\t' -v OFS='\t' '{print $1, $2, $3}' > templates
    templates=$(cat templates)

    echo "$templates" | while read p; do
        templates_query=$(echo "$p" | cut -f1)
        templates_fitness=$(echo "$p" | cut -f2)
        templates_pval=$(echo "$p" | cut -f3)
        helper=$(grep "$templates_query" "$pred_interactome_by_some_pred_bac_dualseq_blast")
        freq=$(echo "$helper" | wc -l)
        templates_fitness2=$(for i in $(seq 1 $freq); do echo "$templates_fitness"; done)
        templates_pval2=$(for i in $(seq 1 $freq); do echo "$templates_pval"; done)
        paste <(echo "$helper") <(echo "$templates_fitness2") <(echo "$templates_pval2");
    done > "$intermediate_files"bacteria_vs_dualseq_values.tsv

fi

if [ ! -f "$pred_interactome_by_some_pred_host_dualseq_blast" ]; then
    blastp -query $pred_interactome_by_some_pred_host  -db filtered_dualseqdb_v1_withseq_human.fasta -evalue 1e-50  -outfmt 6 | awk -F'\t' -v OFS='\t'  '{print $1, $2, $3, $11, $12}' > "$pred_interactome_by_some_pred_host_dualseq_blast"

#human
    query=$(awk -F'\t' -v OFS='\t' '{print $1}' "$pred_interactome_by_some_pred_host_dualseq_blast")
    awk '{print $2}' "$pred_interactome_by_some_pred_host_dualseq_blast" > template
    grep -wf template filtered_dualseqdb_v1_withseq_human_with_scores.tsv | awk -F'\t' -v OFS='\t' '{print $1, $2, $3}' >  templates
    templates=$(cat templates)

    echo "$templates" | while read p; do
        templates_query=$(echo "$p" | cut -f1)
        templates_fitness=$(echo "$p" | cut -f2)
        templates_pval=$(echo "$p" | cut -f3)
        helper=$(grep -w "$templates_query" "$pred_interactome_by_some_pred_host_dualseq_blast")
        freq=$(echo "$helper" | wc -l)
        templates_fitness2=$(for i in $(seq 1 $freq); do echo "$templates_fitness"; done)
        templates_pval2=$(for i in $(seq 1 $freq); do echo "$templates_pval"; done)
        paste <(echo "$helper") <(echo "$templates_fitness2") <(echo "$templates_pval2");
    done > "$intermediate_files"host_vs_dualseq_values.tsv
fi

rm template templates
#host centrality has to be calculated based on the host_dir

if [[ "$host_taxonID" ]]; then

    hi_centrality_name=$(echo hi_centrality_"$host_taxonID".tsv)
    if [[ ! -f "$hi_centrality_name" ]] ; then
        interactome_name_c=$(echo  "$host_taxonID"_interactome.gz)
        proteome_alias_c=$(echo "$host_taxonID"_proteome_alias.gz)
        interactome_name=$(echo "$host_taxonID"_interactome)
        proteome_alias=$(echo "$host_taxonID"_proteome_alias)
        filtered_interactome_name=$(echo filtered_"$host_taxonID"_interactome)
        filtered_proteome_alias=$(echo filtered_"$host_taxonID"_proteome_alias)
        hi_centrality_name=$(echo hi_centrality_"$host_taxonID".tsv)
    # Download host interactome from STRING database with given taxonID
        curl -s "https://stringdb-static.org/download/protein.links.v11.0/$host_taxonID.protein.links.v11.0.txt.gz" > "$interactome_name_c" && gunzip "$interactome_name_c"

# Download alias file

        curl -s "https://stringdb-static.org/download/protein.aliases.v11.0/$host_taxonID.protein.aliases.v11.0.txt.gz" > "$proteome_alias_c" && gunzip "$proteome_alias_c"

# Filter interactome and disccard interactions with combined score < 900

        awk '{ if ($3 >= 900) {print $1, $2} }' "$interactome_name" > "$filtered_interactome_name"


# Parse alias file to get uniprot IDs

        grep 'Ensembl_UniProt_AC' "$proteome_alias" | awk -F '\t' -v OFS='\t' '{print $1, $2}' > "$filtered_proteome_alias"

        time Rscript "$script_dir"/interactome_node_betweenness_centrality.R "$filtered_interactome_name" "$filtered_proteome_alias" "$hi_centrality_name"

    fi
fi


final_interactome_raw_file=$(sed 's/\.tsv$/_raw_with_weighted_scores.tsv/' <<< "$interactome")
final_interactome_with_scores=$(sed 's/\.tsv$/_only_weighted_scores.tsv/' <<< "$interactome")
final_interactome_log2fc_file=$(sed 's/\.tsv$/_help_file_with_logfc.tsv/' <<< "$interactome")
final_interactome_IDs_file=$(sed 's/\.tsv$/_help_file_with_IDS.tsv/' <<< "$interactome")

#echo "$final_interactome_raw_file"
#echo "$final_interactome_with_scores"
#echo "$final_interactome_log2fc_file"
#echo "$final_interactome_IDs_file"



#######

Rscript "$script_dir"/integrate_databases_copy.R "$intermediate_files"host_vs_dualseq_values.tsv "$intermediate_files"bacteria_vs_dualseq_values.tsv "$intermediate_files"bacteria_vs_bacfitbase_values.tsv  "$intermediate_files"bacteria_vs_phi_base_values.tsv "$hi_centrality_name" "$interactome" "$final_interactome_raw_file" "$final_interactome_with_scores" "$final_interactome_log2fc_file" "$final_interactome_IDs_file"

#Rscript "$script_dir"/integrate_databases_copy.R "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_host_dualseq_values.tsv "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_bac_dualseq_values.tsv "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_bac_values.tsv "$host"_vs_"$pathogen"_predicted_interactome_false_positives_"$falsepos_var"_by_"$consensus_predictors"_percent_of_predictors_phi_base_values.tsv "$hi_centrality_name" "$interactome" "$final_interactome_raw_file"  "$final_interactome_with_scores" "$final_interactome_log2fc_file" "$final_interactome_IDs_file"


if [ -f "$custom_desc_list" ]; then
    rm "$custom_desc_list"
fi

