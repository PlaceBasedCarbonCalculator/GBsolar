# install.packages("geosam", repos = c("https://walkerke.r-universe.dev", "https://cloud.r-project.org"))
# usethis::edit_r_environ()

library(geosam)
geosam_install(method = 'conda', conda = "C:/Users/earmmor/AppData/Local/miniconda3/_conda.exe")
geosam_status()
