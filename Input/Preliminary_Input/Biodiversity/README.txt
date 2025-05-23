##################################################
Generation of Input to calculate Species Richness
##################################################


Performance.xlsx and Species_coef.xlsx,i.e., input of the Species Richness calculation, were generated using the following script:

	- InputBiodiveristy1.R # Taxonomic and geographic cleaning of biodiveristy database + thinning (1 records per 10km)
	
	- InputBiodiveristy2.py # Create Pseudo-absence data

	- InputBiodiveristy3.R # Calibrate and validate the SDMs These files can be found in Premilinary_Input/Biodiveristy.


The input of these script can be downloaded at : 10.5281/zenodo.15497375

The input of InputBiodiveristy1.R are :
	- GBIF.xlsx	     data source: GBIF. (2024). Global Biodiversity Information Facility. https://www.gbif.org/
	- SALVE.xlsx         data source: ICMBio. (2024). Sistema de Avaliação do Risco de Extinção da Biodiversidade – SALVE. https://salve.icmbio.gov.br/
	- CRIA.Xlsx          data source: CRIA. (2024). SpeciesLink. https://specieslink.net/

The input of InputBiodiveristy2.R are:
	- HMW_Mammal.gpkg    data source: Marsh, C. J., Sica, Y. V., Burgin, C. J., Dorman, W. A., Anderson, R. C., del Toro Mijares, I., Vigneron, J. G., Barve, V., Dombrowik, V. L., Duong, M., Guralnick, R., Hart, J. A.,  Maypole, J. K., McCall, K., Ranipeta, A., Schuerkmann, A., Torselli, M. A., Lacher, T., Mittermeier, R. A., … Jetz, W. (2022). Expert range maps of global mammal distributions harmonised to three taxonomic authorities. Journal of Biogeography, 49(5), 979–992. https://doi.org/10.1111/jbi.14330

The input of InputBiodiveristy3.R are: 
	- Presence.xlsx
	- Absence.xlsx

