#!/bin/bash

mda_iter_list=(2 3 4)
app_dir_temp="1dthermal_parastate_1month"
sort_obs_inc_list=(0)  # zero or one

for sort_obs_inc in ${sort_obs_inc_list[@]}; do
    for mda_iter in ${mda_iter_list[@]}; do

        # Get the application directory name
        if [ ${sort_obs_inc} -eq 0 ] ; then
            app_dir="${app_dir_temp}_${mda_iter}mda"
        else
            app_dir="${app_dir_temp}_sortobsin_${mda_iter}mda"
        fi
        #app_dir="${app_dir}_mda"
        echo ${app_dir}
        #echo ${assim_reso}
        #echo ${mda_iter}
        #echo ${mda_iter_list[@]}

        python workflow_1dthermal_sh.py ${app_dir} ${mda_iter} ${sort_obs_inc}

        # Check exit code
        if [ $? -ne 0 ]
        then
            echo "Something wrong when running the script..."
            break
        fi

    done
done
