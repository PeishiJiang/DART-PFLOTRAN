#!/bin/bash

#assim_reso_list=(7200  14400)
#assim_reso_list=(18000 21600)
assim_reso_list=(3600)
mda_iter_list=(1)
#mda_iter_list=(1 2 3 4)
#app_dir_temp="1dthermal_test_1month_obserror015_rescaled"
app_dir_temp="1dthermal_test_1month_obserror005"
rescaled=0  # zero or one

for assim_reso in ${assim_reso_list[@]}; do
    for mda_iter in ${mda_iter_list[@]}; do

        # Get the application directory name
        if [ ${rescaled} -eq 0 ] ; then
            app_dir="${app_dir_temp}_${assim_reso}s_${mda_iter}mda"
        else
            app_dir="${app_dir_temp}_rescaled_${assim_reso}s_${mda_iter}mda"
        fi
        #app_dir="${app_dir}_mda"
        echo ${app_dir}
        #echo ${assim_reso}
        #echo ${mda_iter}
        #echo ${mda_iter_list[@]}

        # Run workflow_sh.py
        python workflow_sh.py ${app_dir} ${assim_reso} ${mda_iter} ${rescaled}

        # Check exit code
        if [ $? -ne 0 ]
        then
            echo "Something wrong when running the script..."
            break
        fi

    done
done
