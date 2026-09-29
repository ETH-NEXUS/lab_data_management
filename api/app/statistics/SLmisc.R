
# helper functions
#' A function to concatenate two strings
#'
#' @param a
#' @param b
#'
#' @returns a concatenated string of a and b
#' @export
#' @author Michael Prummer <prummer@nexus.ethz.ch>
#'
#' @examples
'%&%' = function(a,b) {
  paste(a,b,sep="")
}


#' Plot 2D plate heatmaps with ggplot2.
#'
#' @param dat data.frame with columns plate, row, col, and a column with name = type.
#' @param type the name of the data column in dat jsed for plotting.
#' @param plot logical indicating whether the plot should be executed.
#'
#' @returns ggplot object
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
ggplate.hm = function(dat, type = "norm", plot=F){
  # dat - (dataframe) with a factor/character columns plate,
  #       numeric columns row, col
  #       and a numeric column with name = type
  # type - (character) the name of the data column in dat
  # plot - (locigal) indicating whether the plot should be executed

  if(!type %in% c("raw", "norm", "cor", "cpd_type", "logratio")) {
    stop("Plate layout heatmap does not know the type.")
  }
  dat$z = dat[, type]
  if(type != "cpd_type"){
    z.range = quantile(dat$z, c(0.01,0.99), na.rm=T)
    dat$z[dat$z < z.range[1]] = z.range[1]
    dat$z[dat$z > z.range[2]] = z.range[2]
  }
  breaks.x = sort(unique(dat$col))
  breaks.x = breaks.x[seq(1,length(breaks.x),1)]
  range.x = range(dat$col)+c(-0.5,0.5)
  labels.y = LETTERS[sort(unique(dat$row))]
  labels.y = labels.y[seq(1,length(labels.y),1)]
  breaks.y = sort(unique(dat$row))
  breaks.y = breaks.y[seq(1,length(breaks.y),1)]
  range.y = rev(range(dat$row)+c(-0.5,0.5))
  col.pal = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255)

  if(type == "cpd_type"){
    # color palette(s):
    type.pal = brewer.pal(9,"Set1")
    type.pal = c(type.pal[2:1], "black", type.pal[3:9], brewer.pal(8,"Dark2"))
    col.pal = colorRampPalette( rev(brewer.pal(9, "Blues")) )(length(unique(dat$conc)))
    col.pal = c(type.pal[seq_along(unique(dat[, type]))], col.pal)

    p = ggplot(dat, aes(ymin=row-0.5, ymax=row+0.5, xmin=col-0.5, xmax=col+0.5, fill=z)) +
      geom_rect(colour="grey50") + #coord_equal() +
      scale_y_reverse(labels=labels.y, breaks=breaks.y, limits=range.y, expand = c(0,0)) +
      scale_x_continuous(breaks=breaks.x, limits=range.x, expand = c(0,0)) +
      scale_fill_manual(values = col.pal, name="") +
      theme(strip.text = element_text(size=8),
            legend.text=element_text(size=rel(0.7)),
            legend.key.size=unit(3,"mm"),
            axis.text.x = element_text(size=rel(0.7)),
            axis.text.y = element_text(size=rel(0.7))) +
      guides(position="right")
  } else {
    p = ggplot(dat, aes(ymin=row-0.5, ymax=row+0.5, xmin=col-0.5, xmax=col+0.5)) +
      geom_rect(aes(fill=z), colour="grey50", linewidth=0.1) + coord_equal() +
      #geom_tile(aes(fill=z), colour="grey50", size=0.1) + coord_equal() +
      scale_y_reverse(labels=labels.y, breaks=breaks.y, limits=range.y, expand = c(0,0)) +
      scale_x_continuous(breaks=breaks.x, limits=range.x, expand = c(0,0)) +
      #scale_color_brewer(palette="RdYlBu") +
      scale_fill_gradientn(colours=col.pal, name = type) +
      theme(legend.position="right", strip.text = element_text(size=5),
            axis.text.x = element_text(size=rel(0.5)),
            axis.text.y = element_text(size=rel(0.5))) +
      guides(position="top", color=guide_legend(byrow=F, title = NULL))
  }

  if(plot) print(p)
  return(p) # ggplot object
}



#' A function to read data from LDM export.
#'
#' @param path_data path to the data file
#'
#' @returns a data.frame with columns plate, row, col, cpd_type, raw, name, well
#' and more, depending on the content of the LDM export.
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
read.LDM = function(path_data) {
  dat = utils::read.csv(path_data)
  #dat = read.table(path_data, head=T, sep=";")
  # if(ncol(dat)==1) dat = read.table(path_data, head=T, sep=",", as.is = T, comment.char = "")
  if(ncol(dat)==1) dat = read.table(path_data, head=T, sep="\t", comment.char = "")
  if(ncol(dat)==1) stop("Input file has invalid delimiter. Valid format is tsv or csv.")
  id_colnames = c("unique_identifier", "plate", "plate_row", "plate_column", "control", "value")
  #id_colnames = union(id_colnames, names(dat))
  dat = dat[, id_colnames]
  names(dat)[1:6] = c("uid", "plate", "row", "col", "cpd_type", "raw")
  dat$well = sprintf("%s%02.0f", LETTERS[dat$row], dat$col)
  return(dat)
}


#' A function to read data from experiment data.
#'
#' @param path_data path to the data file
#'
#' @returns a data.frame with columns layout, experiment, replicate, cell_type, condition
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
read.ExpData = function(path_data) {
  dat = utils::read.csv(path_data)
  if(ncol(dat)==1) dat = read.table(path_data, head=T, sep=",", as.is = T, comment.char = "")
  if(ncol(dat)==1) dat = read.table(path_data, head=T, sep="\t", comment.char = "")
  if(ncol(dat)==1) stop("Input file has invalid delimiter. Valid format is tsv or csv.")
  if(!"layout" %in% names(dat)) dat$layout = "LAY1"
  if(!"experiment" %in% names(dat)) dat$experiment = "EXP1"
  id_keep = which(dat$plate != "null")
  dat = dat[id_keep, ]
  if(!("replicate" %in% names(dat)) | all(is.na(dat$replicate))) dat$replicate = "RPL1"
  if(!("cell_type" %in% names(dat)) | all(is.na(dat$cell_type))) dat$cell_type = "CTY1"
  if(!("condition" %in% names(dat)) | all(is.na(dat$condition))) dat$condition = "CND1"
  if(!("measurement_label" %in% names(dat)) | all(is.na(dat$measurement_label))) dat$measurement_label = "MEAS01"
  names(dat) = gsub("lib_plate_barcode", "lib_plate", names(dat))
  return(dat)
}



#' A function to read data from a plate layout definition file.
#'
#' @param data.desc.fn path to the data description file
#'
#' @returns a list with the following elements: fdr_cut, act_cut, plate_rows,
#' plate_cols, plate_wells, condi, runs, channels, cpd_desc, meta
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
read.DA_config = function(data.desc.fn) {
  # load data description file
  d.d = readLines(data.desc.fn)
  idx = grep("description", d.d)
  screen.id = idx[2]
  cpd.id = idx[3]
  lay.id = idx[4]
  plt.id = idx[5]
  nrow.dd = length(d.d)
  rm(d.d)
  # load sample description table
  d.head = read.table(data.desc.fn, sep="\t", head=F, skip=screen.id,
                      nrows = cpd.id - screen.id -1, stringsAsFactors = F)
  d.head = d.head[, !apply(d.head, 2, function(x) all(is.na(x)))]
  d.cpd = read.table(data.desc.fn, sep="\t", head=T, skip=cpd.id,
                     nrows = lay.id - cpd.id -1, stringsAsFactors = F)
  d.cpd = d.cpd[,-grep("X", names(d.cpd))]
  d.cpd = d.cpd[-which(is.na(d.cpd$dilution.factor)),]
  d.plt = read.table(data.desc.fn, sep="\t", head=T, skip=plt.id,
                     nrows = nrow.dd - plt.id -1, stringsAsFactors = F)
  d.plt = d.plt[,-grep("X", names(d.plt))]
  #d.plt$plate = as.character(sapply(d.plt$filename, function(x) strsplit(x, "\\.")[[1]][2]))
  d.plt$plate = d.plt$barcode
  d.plt = d.plt[!is.na(d.plt$plate), ]

  l_out = list(fdr_cut = as.numeric(d.head$V3[d.head$V2=="FDR cutoff"]), # FDR cutoff,
               act_cut = as.numeric(d.head$V3[d.head$V2 == "pct-ctr cutoff"]), # fold-change cutoff)
               plate_rows = as.numeric(d.head$V3[d.head$V2=="plate rows"]),
               plate_cols = as.numeric(d.head$V3[d.head$V2=="plate columns"]),
               plate_wells = as.numeric(d.head$V3[d.head$V2=="plate size"]),
               condi = as.character(d.head[which(d.head$V2=="conditions"), -1]),
               runs = as.character(d.head[which(d.head$V2=="replicates"), -1]),
               channels = as.character(d.head[which(d.head$V2=="channels"), -1])
  )
  l_out = lapply(l_out, function(x) x[x!="" & !is.na(x)])
  d.lay = read.table(data.desc.fn, sep="\t", head=F, skip=lay.id,
                     nrows = plt.id - lay.id, stringsAsFactors = F)
  idx = which(apply(d.lay[,-1], 1, function(x) all(x=="")))
  d.lay = d.lay[-idx,-1]
  idx = grep("plate.type", d.lay[,1])
  meta.colnames = d.lay[idx+1,1]
  meta.platype = sapply(idx, function(x) {
    tt = as.character(d.lay[x, -1])
    tt[tt!="" & tt!="NA" & !is.na(tt)]
  } )
  col.numbers = lapply(idx, function(x) as.numeric(d.lay[x+1, -1]))
  t.diff = diff(c(idx,nrow(d.lay)+1))
  row.letters = mapply(function(x,y) d.lay[x:y,1], x=idx+2, y=idx+t.diff-1)
  if(class(row.letters)[1]!="list") {
    row.letters = apply(row.letters, 2, function(x) rbind(as.list(x)))
  }
  row.letters = lapply(row.letters, as.character)
  row.letters = lapply(row.letters, function(x) x[!is.na(x)])
  col.numbers = lapply(col.numbers, function(x) x[!is.na(x)])
  nr.wells.per.plate = sapply(col.numbers,length) * sapply(row.letters, length)
  nr.platypes = length(unique(meta.platype))

  # generate df with correct number of wells for each plate
  meta = data.frame(NULL)
  for(ii in seq(length(meta.colnames))){
    tt = expand.grid(plate.type=meta.platype[ii], row=row.letters[[ii]],
                     col=col.numbers[[ii]], KEEP.OUT.ATTRS=F, stringsAsFactors = F)
    #meta = merge(meta,tt, all.x=T)
    meta = rbind(meta, tt)
  }
  meta = meta[!duplicated(meta),]

  # attach relevant column names from DA_config
  tt = data.frame(name=meta.colnames, plate.type=meta.platype, value=NA, stringsAsFactors=F)
  tt = reshape2::dcast(tt, plate.type~name)
  meta = merge(meta, tt, all.x=T)
  meta$well = sprintf("%s%02.0f", meta$row, meta$col)
  meta$row = match(meta$row, LETTERS)
  meta$uwid = meta$plate %&% meta$well
  # fill variables with values:
  idx = grep("plate.type", d.lay[,1])
  idx.diff = diff(c(idx,nrow(d.lay)+1))
  for(ii in seq(length(idx))){
    #ii=1
    tm = d.lay[(idx[ii]+2):(idx[ii]+idx.diff[[ii]]-1), 2:(length(col.numbers[[ii]])+1)]
    dimnames(tm) = list(row.letters[[ii]], col.numbers[[ii]])
    td = as.data.frame(tm)
    #td$row = rownames(td)
    rd = melt(as.matrix(td), value.name = meta.colnames[[ii]], varnames = c("row", "col"))
    #rd = do.call(rbind, replicate(length(meta.platype[[ii]]), rd, simplify = F))
    rd$plate.type = rep(meta.platype[[ii]], each=nr.wells.per.plate[[ii]])
    rd$uwid = rd$plate.type %&% sprintf("%s%02.0f", rd$row, rd$col)
    meta[match(rd$uwid, meta$uwid), names(rd)[3]] = rd[,3]
  }
  rm(rd)
  # find numeric variables:
  isna = apply(meta,2,function(x) any(is.na(as.numeric(x))))
  meta[, which(!isna)] = apply(meta[, which(!isna), drop=F], 2, as.numeric)
  names(d.plt)[match("type", names(d.plt))] = "plate.type"

  meta = merge(d.plt, meta, all.x=T)

  # str(meta)
  # names(meta)
  # summary(meta)
  meta$type = factor(meta$type, levels=d.cpd$type)
  l_out$cpd_desc = d.cpd
  l_out$meta = meta
  return(l_out)

}



#' A funtion to read library description files.
#'
#' @param path_libfile path to the library description file
#' @param libtype the type of the library description file, one of "actitarg", "nexus_fda",
#' "seleckchem_clinical", "EPC", "LLD"
#'
#' @returns a data.frame with columns CatalogNumber, CompoundName, Target, CAS, MW, and more,
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
read.library.desc = function(path_libfile, libtype = "actitarg"){
  # librype one of "actitarg", "nexus_fda", "seleckchem_clinical", "EPC", "LLD"
  libdesc = "-1"
  if(libtype == "actitarg"){
    libdesc = read.csv(path_libfile)
    libdesc = libdesc[, -c(6:21, 147:149)]
    id_keep = which(!duplicated(libdesc$CatalogNumber))
    libdesc = libdesc[id_keep, ]
  }
  if(libtype == "nexus_fda"){
    libdesc = read.csv(path_libfile)
    libdesc = libdesc[, -c(2, 3, 6:29, 34, 35, 39:52)]
    id_keep = which(!duplicated(libdesc$CatalogNumber))
    libdesc = libdesc[id_keep, ]
  }
  if(libtype == "selleckchem_clinical"){
    libdesc = read.table(path_libfile, sep = "\t", head=T, fill=T,
                         quote = "", comment.char = "")
    libdesc = libdesc[, -c(6, 11:35, 39, 43, 44, 49:52)]
    id_keep = which(!duplicated(libdesc$CatalogNumber))
    libdesc = libdesc[id_keep, ]
    libdesc$CompoundName = gsub("\"", "", libdesc$CompoundName)
    libdesc$Target = gsub("\"", "", libdesc$Target)
  }
  if(libtype == "EPC"){
    libdesc = read.table(path_libfile, sep = "\t", head=T, fill=T,
                         comment.char = "")
    libdesc = libdesc[, -c(6, 12:17, 20:31)]
    id_keep = which(!duplicated(libdesc$CatalogNumber))
    libdesc = libdesc[id_keep, ]
  }
  if(libtype == "LLD"){
    libdesc = read.csv(path_libfile)
    libdesc = libdesc[, -c(2, 6, 12:15, 18:28)]
    id_keep = which(!duplicated(libdesc$CatalogNumber))
    libdesc = libdesc[id_keep, ]
  }

  return(libdesc)
}



#' A helper function to convert a vector to a comma separated string.
#'
#' @param x the vector to be converted
#'
#' @returns a comma separated string
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
vec2cs.string = function(x) {
  paste(x[!is.na(x) & x!=""], collapse=", ")
}



#' A helper function to concatenate strings to a maximum length.
#'
#' @param chr_vec the vector of strings to be concatenated
#' @param nchar_max the maximum number of characters
#'
#' @returns a concatenated string with a maximum number of characters
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
concat.chr.max = function(chr_vec, nchar_max = 80) {
  iinchar = nchar(chr_vec)
  iinchar = cumsum(iinchar + 2)
  tmp = chr_vec[which(iinchar <= nchar_max)]
  if(!identical(tmp, chr_vec)) {
    chr_vec = c(tmp, "...")
  }
  return(vec2cs.string(chr_vec))
}



#' A function to generate a table with column names and possible values.
#'
#' @param meta a data.frame with columns and values
#' @param nchar_max the maximum number of characters
#'
#' @returns a table with column names and possible values with a maximum number of characters
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
colNames.colValues = function(meta, nchar_max = 80) {
  tab = matrix(NA, nrow=ncol(meta), ncol=2)
  colnames(tab) = c("column name", "possible values")
  tab[,1] = names(meta)
  l_vals = apply(meta,2,function(x) sort(unique(x)))
  for(ii in seq_along(l_vals)){
    iinchar = nchar(l_vals[[ii]])
    iinchar = cumsum(iinchar + 2)
    tmp = l_vals[[ii]][which(iinchar <= nchar_max)]
    if(!identical(tmp, l_vals[[ii]])) {
      l_vals[[ii]] = c(tmp, "...")
    }
  }
  l_vals
  tab[,2] = as.character(sapply(l_vals, vec2cs.string))
  return(tab)
}


#' A function to compute median and MAD for each plate and cpd_type.
#'
#' @param dd data.frame with columns plate, cpd_type, meas, and a column with name = signal
#' @param signal the name of the data column in dd used for plotting.
#'
#' @returns a data.frame with columns plate, cpd_type, meas, med, and mad
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
get_plateStats = function(dd, signal){
  plate.median = tapply(dd[, signal], list(dd$plate, dd$cpd_type, dd$meas), median, na.rm=T)
  plate.median = melt(plate.median)
  names(plate.median) = c("plate", "cpd_type", "meas", "med")
  plate.mad = tapply(dd[, signal], list(dd$plate, dd$cpd_type, dd$meas), mad, na.rm=T)
  plate.mad = melt(plate.mad)
  names(plate.mad) = c("plate", "cpd_type", "meas", "mad")
  plate_stats = inner_join(plate.median, plate.mad)

  return(plate_stats)
}



#' A function to fit a 4 parameter Hill equation model to a set of dose response data.
#'
#' @param data a data.frame in long format, with columns name, conc, response
#'
#' @returns a list of drm fit objects
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
fit.DRC = function(data){
  subset_data = split(data, data$name)
  y = lapply(subset_data, function(x) {
    x$logconc = log10(x$conc)
    drm(response ~ conc, data = x, fct = LL2.4())
  })
  names(y) = names(subset_data)
  return(y)
}


#' A function to extract the parameters of a drm fit object.
#'
#' @param model_fits a list of drm fit objects from 'fit.DRC()'
#'
#' @returns a list of data.frames with the parameters of the drm fit objects
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
get.fit.params = function(model_fits){
  y <- lapply(model_fits, function(model) {
    out = summary(model)$coef
    out = cbind(out, confint(model))
    tmp = gsub("(.*):.Intercept.", "\\1", rownames(out))
    rownames(out) = tmp
    out = cbind(term=tmp, as.data.frame(out))
    return(out)
  })
  return(y)
}


#' Return a dataframe with fit parameters from a dose response curve
#'
#' @param model_params a list of data.frames with the parameters of the drm fit objects, from 'get.fit.params()'
#'
#' @returns a data.frame with columns name, EC50, EC50_se, cilo, cihi, plalo, plahi, plalo_se, plahi_se, halfamp
#' @export
#' @author Michael Prummer michael.prummer@nexus.ethz.ch
#'
#' @examples
tidy.DRC = function(model_params){
  EC50 = sapply(model_params, function(x) x$Estimate[x$term=="e"])
  cilo = sapply(model_params, function(x) x$'2.5 %'[x$term=="e"])
  cihi = sapply(model_params, function(x) x$'97.5 %'[x$term=="e"])
  stderr = sapply(model_params, function(x) x$'Std. Error'[x$term=="e"])
  plalo = sapply(model_params, function(x) x$Estimate[x$term=="c"])
  plahi = sapply(model_params, function(x) x$Estimate[x$term=="d"])
  # Extract the std-err of EC50, lower plateau, and upper plateau values for each compound
  EC50_stderr <- sapply(model_params, function(x) x$'Std. Error'[x$term == "e"])
  lower_plateau_stderr <- sapply(model_params, function(x) x$'Std. Error'[x$term == "c"])
  upper_plateau_stderr <- sapply(model_params, function(x) x$'Std. Error'[x$term == "d"])
  d_ec50 = data.frame(name = names(model_params), EC50=exp(EC50), EC50_se=exp(EC50_stderr),
                      cilo=exp(cilo), cihi=exp(cihi), plalo=plalo, plahi=plahi,
                      plalo_se=lower_plateau_stderr, plahi_se=upper_plateau_stderr,
                      halfamp = (plalo + plahi)/2)
  return(d_ec50)
}



# 1. Separate plate heatmap panel for each channel raw signal
# 2. Channel trafo options: none (y=x), log (y=log10(x))
# 3. Channel assignemnt: y1, y2
# 4. Channel combination options: none (z=y1), ratio (z=y1/y2), difference (z=y1-y2)
# 5. Normalization options: none (norm=z), Z-score (norm=(z-median(z_neg))/sd(z_neg)),
#                           pct-ctr (norm=(z-median(z_neg))/(median(z_pos) - median(z_neg)))
# 6. Correction options: none, polynomial, median polish
# 7. Multiple measurements - conditions: none, groups, series
# 8. Multiple measurements - replicates: Number (0, 1, 2, ...)

# ccc_ip - cell count correction factor per well in plate:
#   ccc_ip = cc_ip/<cc_p>
# <cc_p> - geometric average of cell count signal over all wells of a plate

# z_norm_ip = log(z_ip) - log(z_neg_p) - log(ccc_ip)
# z_norm_ip ~ condi_p = Z-test: (z_norm_ip(C1) - z_norm_ip(C2)) / sqrt(sd_ip(C1)*sd_ip(C2)) == 0 ?


