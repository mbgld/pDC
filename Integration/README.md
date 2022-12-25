# Integration of 10 datasets comprised of tumor and normal tissues and annotation for comparision of TME of current smoker's lung adenocarcinoma with never smoker's.

Combining indivisual dataset 
Integrate by canonical correlation anaysis in seurat package and clustering.
Annotate each clusters. For cancer cell annotation use the barcodes of cancer cell obtained from 'each_case' folder.  For other cells, annotate them comparing FindMarker's results and known cannonical markers.
Split the integrate seurat object into major cell clusters for additional in depth analysis.
