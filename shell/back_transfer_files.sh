rsync -av -e ssh --exclude-from 'shell/exclude_list.txt' lbenz@fasselogin.rc.fas.harvard.edu:~/tv_effects/* .
#Rscript scripts/util/img_to_pdf.R
#rm Rplots.pdf