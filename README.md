# MRCode
twosample MR code
suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(gwasvcf)
    library(gwasglue)
    library(ieugwasr)
    library(dplyr)
    library(curl)
    library(parallel)
    library(rvest)  # 爬 gwas info
    library(foreach)
    library(future)
    library(future.apply)
})

options(future.globals.maxSize = 20 * 1024^3)

opt <- list()
opt$exposure_id_list <- "Non-ischemic cardiomyopathy:finn-b-I9_NONISCHCARDMYOP"
opt$outcome_id_list <- "Chronic kidney disease:ebi-a-GCST003374"
opt$vcf_download_method <- "manu"
opt$p_value =  5e-08
opt$kb = 10000
opt$r2 = 0.001
opt$vcf_dir_name = "VCF"
opt$output_dir = "06_Reverse_MR_new"
opt$token = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxMjY5ODA3Njc2QHFxLmNvbSIsImlhdCI6MTczNDgzNzc5OSwiZXhwIjoxNzM2MDQ3Mzk5fQ.HI1VWKlOdnUe9dd0hM58uBrji0kxNHOVnEFffbUSftJHd3vsg8TQyxkLR_89tvafcuqHEVPJm3u1kH-V-yt_YtYEBSEeH3KWX6k7mS23gQsHJtatnrcWKKZ_gL8ZXL2HKUH_9mewa-em9_8jsNfTugYPlq2jkNlgitoXgJaswtkOdMuH1WtvM9q9Fm5FywzcVNDH5GzqJJhqdp4c-8eNuxmCFr4UAoTQTg9iQ63td-Sn8hON14XSsxDjTO42Z2HCnWjZrDMel5PxMx_nK9KHnQ278vNVkQuJDNGfo9mw4RxOxAnq0tHA98TW7EdDwgm3z5_WU_szRZ2W2C2zOVYg_w"


red <- "\033[31m"
green <- "\033[32m"
yellow <- "\033[33m"
blue <- "\033[34m"
magenta <- "\033[35m"
cyan <- "\033[36m"
reset <- "\033[0m"

# colors
colors <- c('#4DBBD5FF','#E64B35FF','#00A087FF','#3C5488FF','#925E9FFF','#91D1C2FF',
            '#8491B4FF','#7E6148FF','#0072B5FF','#E18727FF','#B09C85FF','#20854EFF',
            '#6F99ADFF','#FFDC91FF','#00468BFF','#FDAF91FF','#B24745FF','#6699FFFF',
            '#99991EFF','#FFCCCCFF','#358000FF','#99CCFFFF','#FFCC00FF','#8C564BFF',
            '#BCBD22FF','#996600FF','#5CB85CFF','#F39B7FFF','#CE3D32FF','#749B58FF',
            '#466983FF','#F0E685FF','#D595A7FF','#924822FF','#7A65A5FF','#C75127FF',
            '#FFA319FF','#8A9045FF','#8F3931FF','#00AF66FF','#748AA6FF','#D0DFE6FF',
            '#C71000FF','#008EA0FF','#8A4198FF','#D5E4A2FF','#5A9599FF','#FF6348FF',
            '#B7E4F9FF','#FF95A8FF','#526E2DFF','#FB6467FF','#E89242FF','#69C8ECFF',
            '#917C5DFF','#FED439FF','#709AE1FF','#D2AF81FF','#FD7446FF','#4DBBD5FF',
            '#E64B35FF','#00A087FF','#3C5488FF','#925E9FFF','#91D1C2FF','#8491B4FF')


get_gwas_info <- function(gwas_id) {
  
  url <- paste0("https://gwas.mrcieu.ac.uk/datasets/", gwas_id, "/")
  
  
  page <- tryCatch({
    read_html(url)
  }, error = function(e) {
    message("Error reading URL: ", url)
    return(NULL)
  })
  

  if (is.null(page)) return(NULL)
  
  
  table <- page %>% html_node("table")
  
  if (is.null(table)) {
    message("No table found for ", gwas_id)
    return(NULL)
  }
  
  
  rows <- table %>% html_nodes("tr")
  gwas_info <- sapply(rows, function(row) {
    key <- row %>% html_node("th") %>% html_text(trim = TRUE) 
    value <- row %>% html_node("td") %>% html_text(trim = TRUE)  
    if (!is.null(key) && !is.null(value)) {
      return(c(key, value))  # 返回键值对
    } else {
      return(NULL)
    }
  })
    
  gwas_info <- as.data.frame(gwas_info, stringsAsFactors = FALSE)  
  gwas_info <- setNames(gwas_info[-1, ], gwas_info[1, ])  

  
  cols_to_process <- c("Year", "ncase", "ncontrol", "Sample size", "Number of SNPs")
  select_cols <- colnames(gwas_info)[colnames(gwas_info) %in% cols_to_process]
  gwas_info[select_cols] <- lapply(gwas_info[select_cols], function(x) gsub(",", "", x))
  gwas_info <- gwas_info %>% mutate(across(all_of(select_cols), as.numeric))
  is_valid <- !any(grepl("Sample size", colnames(gwas_info))) &&  any(grepl("ncase", colnames(gwas_info))) && any(grepl("ncontrol", colnames(gwas_info)))
  if( is_valid ){gwas_info[["Sample size"]] <- gwas_info$ncase + gwas_info$ncontrol}
  
  return(gwas_info)
}

id_get_download_path <- function(ids = c(exposure_id_list_value, outcome_id_list_value), finngen_manifest_path = "/mnt/f/Data/Biodocument/MR/finngen_R12/finngen_R12_manifest.tsv" ){
        id_download_paths <- lapply(ids,function(id){
                            if(grepl("^finn",id)){
                              
                                
                                http_data <- read.delim(finngen_manifest_path)
                                id_name <- c(exposure_id_list_name,outcome_id_list_name)[which(c(exposure_id_list_value,outcome_id_list_value) %in% id)]
                                id_phenotypes <- http_data$phenotype[grep(tolower(id_name),tolower(http_data$phenotype))]
                                id_http_paths <- http_data$path_https[grep(tolower(id_name),tolower(http_data$phenotype))]
                                metadata <- data.frame(phenotype = id_phenotypes,http_path = id_http_paths)
                                if( length(id_phenotypes) == 0 ){
                                   
                                   if(!grepl(" ",id_name)){ 
                                       cat("\n",yellow,paste0("Phenotype: ",blue,id_name),reset,"\n")
                                       cat("\n",yellow,"Please input a special character: ",reset,"\n")
                                       flush.console()  
                                       special_character <- readline(prompt = "")
                                       id_name <- gsub(special_character," ",id_name)
                                   }
                                    
                                   
                                   id_name1 <- id_name %>% strsplit(" ") %>% unlist() %>% .[-length(.)] %>% paste0(.,collapse = " ")
                                   id_name2 <- id_name %>% strsplit(" ") %>% unlist() %>% .[-1] %>% paste0(.,collapse = " ")
                                   id_names <- paste(c(id_name1,id_name2),collapse = "|")
                                   id_phenotypes <- http_data$phenotype[grep(tolower(id_names),tolower(http_data$phenotype))]
                                   id_http_paths <- http_data$path_https[grep(tolower(id_names),tolower(http_data$phenotype))]
                                   metadata <- data.frame(phenotype = id_phenotypes,http_path = id_http_paths)
                                   if(length(id_phenotypes) == 0){
                                       
                                       id_names <- id_name %>% strsplit(" ") %>% unlist() %>% paste0(.,collapse = "|")
                                       id_phenotypes <- http_data$phenotype[grep(tolower(id_names),tolower(http_data$phenotype))]
                                       id_http_paths <- http_data$path_https[grep(tolower(id_names),tolower(http_data$phenotype))]
                                       metadata <- data.frame(phenotype = id_phenotypes,http_path = id_http_paths)
                                       if(length(id_phenotypes) == 1){
                                            id_http_path <- id_http_paths
                                       }else if(length(id_phenotypes) == 0){
                                            stop(paste0("\n",yellow,"Please confirm the : ",blue,id_name,yellow,"    No phenotype matched in FinnGen R12.",reset,"\n"))
                                       }else{
                                            
                                          
                                            cat("\n",yellow, "Column:", "phenotype", reset,"\n")
                                            for (phenotype in metadata[['phenotype']]) {  print(phenotype) }
                                            cat("\n",yellow,paste0("Phenotype: ",blue,id_name),reset,"\n")
                                            cat("\n",yellow,"Please select phenotype:",reset,"\n")
                                            flush.console() 
                                            
                                            phenotype_select <- readline(prompt = "")
                                            id_http_path <- metadata$http_path[metadata$phenotype == phenotype_select]
                                       }
                                       
                                   }else if(length(id_phenotypes) == 1){
                                       
                                       id_http_path <- id_http_paths
                                   }else{
                                        
                                        
                                        cat("\n",yellow, "Column:", "phenotype", reset,"\n")
                                        for (phenotype in metadata[['phenotype']]) {  print(phenotype) }
                                        cat("\n",yellow,paste0("Phenotype: ",blue,id_name),reset,"\n")
                                        cat("\n",yellow,"Please select phenotype:",reset,"\n")
                                        flush.console()  
                                        
                                        phenotype_select <- readline(prompt = "")
                                        id_http_path <- metadata$http_path[metadata$phenotype == phenotype_select]
                                   }
                        
                                }else if(length(id_phenotypes) == 1){
                                    
                                    id_http_path <- id_http_paths
                                }else{
                                    
                                    
                                    cat("\n",yellow, "Column:", "phenotype", reset,"\n")
                                    for (phenotype in metadata[['phenotype']]) { print(phenotype) }
                                    cat("\n",yellow,paste0("Phenotype: ",blue,id_name),reset,"\n")
                                    cat("\n",yellow,"Please select phenotype:",reset,"\n")
                                    flush.console()  
                                    
                                    phenotype_select <- readline(prompt = "")
                                    id_http_path <- metadata$http_path[metadata$phenotype == phenotype_select]
                                }
                            }else{
                                
                                id_http_path <- file.path("https://gwas.mrcieu.ac.uk/files",id,paste0(id,".vcf.gz"))
                            }
                        
                            if(grepl(id,exposure_id_list_value) %>% any()){ id_type <- "exposure" }else{ id_type <- "outcome" }
                            cat("\n",yellow,paste0(id_type,": ",green,id,reset,"\n\t",id_http_path,"\n"))
                            flush.console()  
                        
                            return(id_http_path)
                        })
    id_download_paths <- setNames(id_download_paths, c(exposure_id_list_value,outcome_id_list_value))
    
    return(id_download_paths)
}


vcf_data_download <- function(download_method = opt$vcf_download_method, id_paths = id_download_paths, vcf_save_path = file.path(opt$output_dir,opt$vcf_dir_name) ){
    
    download_file <- function(url, save_path) { system(paste("wget -O", save_path, url)) }
    save_paths <- file.path(vcf_save_path,paste0( names(id_paths), ".vcf.gz"))
    
    if(tolower(download_method) == "auto"){
      
        cat("\n",yellow,"Data Downloading, please waiting...",reset,"\n")
        flush.console()  
        
        
        mclapply(1:length(id_paths), function(i) {
              download_file(id_paths[[i]], save_paths[i])
        }, mc.cores = 15)
        
        cat("\n",yellow,"Data Download End!",reset,"\n")
        flush.console()  
                
    }else{
        
        cat("\n",yellow,"Please manually download the VCF data using the provided ID and download URL, and rename it to '<id>.vcf.gc' format.",reset,"\n")
        cat("\n",yellow,"Path where the VCF file is saved: ",reset,vcf_save_path,"\n")
        cat("\n",yellow,"Have all the exposure and outcome data been downloaded?",reset,"\n")
        flush.console() 
        num = 1
        while(num == 1){
            download_finish <- readline(prompt = "yes?")
            if(download_finish == "yes"){
                
                download_finish_vcf_files <- list.files(vcf_save_path,pattern = ".vcf.gz",full.names = T)
                if(all(download_finish_vcf_files %in% save_paths)){ num = 0 }
            }
        }
    }
}



read_data_and_col_rename <- function(id, file_path = file.path(opt$output_dir, opt$vcf_dir_name)){
    vcf_file_path <- file.path(file_path,paste0(id,".vcf.gz"))
    ncase <- data[[id]]$gwasinfo$ncase
    ncontrol <- data[[id]]$gwasinfo$ncontrol
    sample_size <- data[[id]]$gwasinfo$`Sample size`
    
    
    if( grepl("^finn",id) ){
        vcf_data <- vroom::vroom(vcf_file_path, col_names = T)
        mr_data <- vcf_data %>% rename('chr.exposure' = '#chrom',
                                       'pos.exposure' = 'pos',
                                       'other_allele.exposure' = 'ref',
                                       'effect_allele.exposure' = 'alt',
                                       'SNP' = 'rsids',
                                       'pval.exposure' = 'pval',
                                       'beta.exposure' = 'beta',
                                       'se.exposure' = 'sebeta',
                                       'eaf.exposure' = 'af_alt'
                                       ) %>%
                select('chr.exposure','pos.exposure','effect_allele.exposure', 'other_allele.exposure',
                           'SNP','pval.exposure','eaf.exposure','beta.exposure','se.exposure') %>%
                mutate(pval.exposure = as.numeric(pval.exposure),
                       exposure = id,
                       id.exposure = id,
                       #id.exposure = outcome_id,  
                       samplesize.exposure = sample_size,  
                       ncase.exposure = ncase, 
                       ncontrol.exposure = ncontrol 
                      )
    }else{
        vcf_data <- VariantAnnotation::readVcf(vcf_file_path)
        mr_data <- gwasglue::gwasvcf_to_TwoSampleMR(vcf_data)
    }

    return(mr_data)
}

# CMplot
CMplot_function <- function(mr_data = mr_data, output_dir, p_value = opt$p_value, num_cores = 15) {
    
    two_sample_MR_data <- mr_data[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
    colnames(two_sample_MR_data) <- c("SNP", "CHR", "BP", "pvalue")

    current_path <- getwd()
    setwd(output_dir)

   
    plan(multicore, workers = num_cores)

    plot_types <- c("c", "m")
    plot_formats <- c("tiff", "pdf")
    plot_combinations <- expand.grid(type = plot_types, format = plot_formats)

   
    future_lapply(1:nrow(plot_combinations), function(i) {
        type <- as.character(plot_combinations$type[i])
        format <- as.character(plot_combinations$format[i])

        if (type == "m") width <- 8 else width <- 6
        CMplot::CMplot(two_sample_MR_data,
                       plot.type = type,
                       LOG10 = TRUE,
                       threshold = p_value,
                       threshold.lty = 1, threshold.lwd = 1.5,
                       cex = 0.2, signal.cex = 0.4,
                       chr.den.col = NULL,
                       bin.size = 1e5,
                       file = format,
                       dpi = 300,
                       height = 6,
                       width = width,
                       ylim = c(0, 60))
    })

    plan(sequential)
    
    setwd(current_path)
}


mr_data_clump <- function(mr_data, population, kb = opt$kb, r2= opt$r2, p_value = opt$p_value,
                          bfile_path = "/mnt/f/Data/Biodocument/MR/bfile",
                          plink_bin = "/mnt/f/Data/Biodocument/MR/plink_bin/plink" ){
  
    mr_data_filter <- subset(mr_data, pval.exposure < p_value)

 
    if(population == "European"){bfile_end = "EUR"}
    if(population == "African"){bfile_end = "AFR"}
    if(population == "Admixed American"){bfile_end = "AMR"}
    if(population == "East Asian"){bfile_end = "EAS"}
    if(population == "South Asian"){bfile_end = "SAS"}
    if(population == "Mixed"){bfile_end = "MIX"}
    if(population == "Hispanic or Latin American"){bfile_end = "AMR"}
    

    
    mr_data_clumped <- ld_clump(
         
          dat = tibble(rsid = mr_data_filter$SNP,  
                       pval = mr_data_filter$pval.exposure,  
                       id = mr_data_filter$exposure), 
          clump_kb = kb,   
          clump_r2 = r2,     
          clump_p = 1,         
          bfile = file.path(bfile_path,bfile_end),    
          plink_bin =  plink_bin 
    )
    
    mr_data_clumped <- mr_data_filter %>% filter(mr_data_filter$SNP %in% mr_data_clumped$rsid)
    
    return(mr_data_clumped)
}

# MR 分析出图
MR_Ananysis_Plot <- function(data = df, output_dir = file.path(opt$output_dir,outcome_id,exposure_id)){
  
    mr_result <- mr(data)
    write.csv(mr_result %>% select(-c("id.exposure","id.outcome")),file = file.path(output_dir,"01_mr_result.csv"),row.names = F)

   
    result_or <- generate_odds_ratios(mr_result %>% select(-c("id.exposure","id.outcome")))
    write.csv(result_or,file = file.path(output_dir,"01_mr_result_or.csv"),row.names = F)
    

    heterogeneity_result <- mr_heterogeneity(data)
    write.csv(heterogeneity_result %>% select(-c("id.exposure","id.outcome")),file = file.path(output_dir,"02_heterogeneity_result.csv"),row.names = F)

    
    pleiotropy_result <- mr_pleiotropy_test(data)
    write.csv(pleiotropy_result %>% select(-c("id.exposure","id.outcome")),file = file.path(output_dir,"03_pleiotropy_result.csv"),row.names = F)

    
    directionality_test_result <- directionality_test(data)
    write.csv(directionality_test_result %>% select(-c("id.exposure","id.outcome")),file = file.path(output_dir,"04_directionality_test_result.csv"),row.names = F)
    
    steiger_sl <- steiger_filtering(data)
    write.csv(steiger_sl %>% select(-c("id.exposure","id.outcome")),file = file.path(output_dir,"04_steiger_filtering_result.csv"),row.names = F)

    
    p1 <- mr_scatter_plot(mr_result, data)
    tiff(filename = file.path(output_dir,"05_scatter_plot.tiff"),width = 6,height = 6.5,units = "in",compression = "lzw",res = 400,type = "cairo")
    print(p1)
    dev.off()
    
    pdf(file = file.path(output_dir,"05_scatter_plot.pdf"),width = 6,height = 6.5)
    print(p1)
    dev.off()

   # 森林图
    result_single <- mr_singlesnp(data)
    p2 <- mr_forest_plot(result_single)
    tiff(filename = file.path(output_dir,"06_forest_plot.tiff"),width = 7,height = 6 + nrow(result_single)*0.05,units = "in",compression = "lzw",res = 400,type = "cairo")
    print(p2)
    dev.off()
    
    pdf(file = file.path(output_dir,"06_forest_plot.pdf"),width = 7,height = 6 + nrow(result_single)*0.05)
    print(p2)
    dev.off() 

    # 留一图
    result_loo <- mr_leaveoneout(data)
    p3 <- mr_leaveoneout_plot(result_loo)
    tiff(filename = file.path(output_dir,"07_leaveoneout_plot.tiff"),width = 6.5,height = 6 + nrow(result_loo)*0.05,units = "in",compression = "lzw",res = 400,type = "cairo")
    print(p3)
    dev.off()
    
    pdf(file = file.path(output_dir,"07_leaveoneout_plot.pdf"),width = 6.5,height = 6 + nrow(result_loo)*0.05)
    print(p3)
    dev.off()

    # 漏斗图
    result_single <- mr_singlesnp(data)
    p4 <- mr_funnel_plot(result_single)
    tiff(filename = file.path(output_dir,"08_funnel_plot.tiff"),width = 6,height = 6,units = "in",compression = "lzw",res = 400,type = "cairo")
    print(p4)
    dev.off()
    
    pdf(file = file.path(output_dir,"08_funnel_plot.pdf"),width = 6,height = 6)
    print(p4)
    dev.off()
}

# 解析参数
exposure_id_list <- opt$exposure_id_list %>% strsplit(";") %>% unlist() %>% strsplit(":")  %>% sapply(function(x) x[2]) %>% strsplit(",")
exposure_id_list_value <- exposure_id_list %>% unlist()                                                                                             
names(exposure_id_list) <- opt$exposure_id_list %>% strsplit(";") %>% unlist() %>% strsplit(":")  %>% sapply(function(x) x[1])
exposure_id_list_name <- rep(names(exposure_id_list), times = exposure_id_list %>% sapply(length))                                                                                                       

outcome_id_list <- opt$outcome_id_list %>% strsplit(";") %>% unlist() %>% strsplit(":")  %>% sapply(function(x) x[2]) %>% strsplit(",")
outcome_id_list_value <- outcome_id_list %>% unlist()
names(outcome_id_list) <- opt$outcome_id_list %>% strsplit(";") %>% unlist() %>% strsplit(":")  %>% sapply(function(x) x[1])
outcome_id_list_name <- rep(names(outcome_id_list), times = outcome_id_list %>% sapply(length))                                                                                                                                                                                                                  

opt$exposure_id_list <- exposure_id_list
opt$outcome_id_list <- outcome_id_list

str(opt)

# 创建保存VCF及分析结果的文件夹
dir_paths <- file.path(opt$output_dir,outcome_id_list_value,exposure_id_list_value)
lapply(dir_paths,function(dir_path){
    if(!dir.exists(dir_path)){dir.create(dir_path,recursive = T)}
})
if(!dir.exists(file.path(opt$output_dir,opt$vcf_dir_name))){dir.create(file.path(opt$output_dir,opt$vcf_dir_name),recursive = T)}

# 获取暴露和结局的下载地址
id_download_paths <- id_get_download_path(ids = c(exposure_id_list_value, outcome_id_list_value),
                                          finngen_manifest_path = "/home/zuoanjian/bioinfo_documents/MR/finngen_R12/finngen_R12_manifest.tsv" )

# 创建data对象并获取 gwas info 
data <- vector("list", length = length(id_download_paths))
names(data) <- c(exposure_id_list_value,outcome_id_list_value)

for(name in names(data) ){
    data[[name]][['http_path']] <- id_download_paths[[name]]  # 添加 http_path
    data[[name]][['gwasinfo']] <- get_gwas_info(gwas_id = name)  # 爬虫获取 gwas info 
    #data[[name]][['gwasinfo']] <- ieugwasr::gwasinfo(id = name)  # 需要使用令牌
    if(name %in% exposure_id_list_value){
        data[[name]][['type']] <- "exposure"
    }else if(name %in% outcome_id_list_value){
        data[[name]][['type']] <- "outcome"
    }
}

# 下载 vcf 数据 
#vcf_data_download()

plan(multicore, workers = 1)

exposure_ids <- exposure_id_list_value
outcome_ids <- outcome_id_list_value
ids_combinations <- expand.grid(outcome_ids = outcome_ids, exposure_ids = exposure_ids)

future_lapply(1:nrow(ids_combinations), function(i) {
    outcome_id <- as.character(ids_combinations$outcome_ids[i])
    exposure_id <- as.character(ids_combinations$exposure_ids[i])
    
    cat("\n",yellow,"start: ",exposure_id," to ",outcome_id,reset,"\n")
    flush.console()  # 强制刷新控制台输出


    exposure_data <- read_data_and_col_rename(id = exposure_id)
    outcome_data <- read_data_and_col_rename(id = outcome_id)

    CMplot_function(mr_data = exposure_data, output_dir = file.path(opt$output_dir,outcome_id,exposure_id), p_value = opt$p_value)


    exposure_data_clumped <- mr_data_clump(mr_data = exposure_data, population = data[[exposure_id]]$gwasinfo$Population,
                                           bfile_path = "/home/zuoanjian/bioinfo_documents/MR/bfile",
                                           plink_bin = "/home/zuoanjian/bioinfo_documents/MR/plink_bin/plink"
                                          )
    write.csv(exposure_data_clumped,file = file.path(opt$output_dir,outcome_id,exposure_id,paste0("00_exposure_data_SNP",dim(exposure_data_clumped)[1],".csv")))

  
    data_common <- outcome_data[outcome_data$SNP %in% exposure_data_clumped$SNP,]


    outcome_df <- format_data(outcome_data,
                              type = "outcome", 
                              snps = data_common$SNP, 
                              snp_col = "SNP",  
                              beta_col = "beta.exposure",  
                              se_col = "se.exposure", 
                              eaf_col = "eaf.exposure",
                              effect_allele_col = "effect_allele.exposure",
                              other_allele_col = "other_allele.exposure", 
                              pval_col = "pval.exposure",
                              samplesize_col = "samplesize.exposure",
                              ncase_col = "ncase.exposure", 
                              id_col = "exposure")
    outcome_df$outcome <- outcome_df$id.outcome
    write.csv(outcome_df,file = file.path(opt$output_dir,outcome_id,exposure_id,paste0("00_outcome_data_SNP",dim(outcome_df)[1],".csv")))
    
 
    df <- harmonise_data(exposure_dat = exposure_data_clumped, outcome_dat = outcome_df)
    write.csv(df,file = file.path(opt$output_dir,outcome_id,exposure_id,"00_harmonise_data.csv"),row.names = F)

    data[[exposure_id]][["harmonise_data"]] <- df

    # MR 分析出图
    MR_Ananysis_Plot(data = df, output_dir = file.path(opt$output_dir,outcome_id,exposure_id))
},future.seed = TRUE)

plan(sequential)

saveRDS(data,file = file.path(opt$output_dir,outcome_id,paste0("exposure",length(exposure_id_list_value),"_outcome",length(outcome_id_list_value),".rds")))
