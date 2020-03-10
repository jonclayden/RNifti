#! /usr/bin/env Rscript

library(ore)

values <- symbols <- argPieces <- NULL

for (file in file.path("inst","include","niftilib",c("nifti1_io.h","nifti2_io.h"))) {
    lines <- readLines(file)
    
    wrapped <- which(lines %~% "^[^\\(]+\\)")
    lines[wrapped-1] <- paste(lines[wrapped-1], lines[wrapped])
    lines <- lines[-wrapped]
    
    pieces <- groups(ore.search("^(\\w[\\w\\s]+?\\**)\\s*(\\w+)\\s*\\(([^\\)]+)\\)\\s*;\\s*$", lines))
    values <- c(values, ore.subst("^\\s+(.+?)\\s+", "\\1", pieces[,1,]))
    symbols <- c(symbols, pieces[,2,])
    
    args <- ore.subst("^\\s*(.+?)\\s*$", "\\1", pieces[,3,])
    args <- ore.subst("\\s+([\\s,])", "\\1", args, all=TRUE)
    argPieces <- c(argPieces, groups(ore.search("(\\Avoid\\Z)|\\b([\\w\\s]+?\\**)\\s*(\\w+)\\s*(\\[\\])?(,|\\Z)", args, all=TRUE)))
}

values <- ore.subst("\\s+(\\s)", "\\1", values, all=TRUE)

duplicated <- duplicated(symbols)
abbreviatedSymbols <- ore.subst("nifti", "nii", symbols)
voidType <- sapply(argPieces, function(x) isTRUE(!is.na(x[,1])))
argTypes <- ifelse(voidType, "void", NA_character_)
argTypes[!voidType] <- sapply(argPieces[!voidType], function(x) paste(x[,2],collapse=", "))
argNames <- sapply(argPieces, function(x) paste(x[,3],collapse=", "))
args <- sapply(argPieces, function(x) paste(x[,2],x[,3],collapse=", "))

registerCallables <- paste0("R_RegisterCCallable(\"RNifti\", \"", abbreviatedSymbols[!duplicated], "\", (DL_FUNC) &", symbols[!duplicated], ");")
declarePointers <- paste0("static ", values[!duplicated], "(*_", symbols[!duplicated], ")(", argTypes[!duplicated], ") = NULL;")
mapPointers <- paste0("_", symbols[!duplicated], " = (", values[!duplicated], "(*)(", argTypes[!duplicated], ")) R_GetCCallable(\"RNifti\", \"", abbreviatedSymbols[!duplicated], "\");")

defineWrappers <- rep(NA, length(argPieces))
defineWrappers[voidType & !duplicated] <- paste0(values[voidType & !duplicated], " ", symbols[voidType & !duplicated], " () { NIFTILIB_WRAPPER_BODY_VOID(_", symbols[voidType & !duplicated], ") }")
defineWrappers[!voidType & !duplicated] <- paste0(values[!voidType & !duplicated], " ", symbols[!voidType & !duplicated], " (", args[!voidType & !duplicated], ") { NIFTILIB_WRAPPER_BODY(_", symbols[!voidType & !duplicated], ", ", argNames[!voidType & !duplicated], ") }")
defineWrappers <- na.omit(defineWrappers)

replacePlaceholder <- function (file, labels, strings) {
    lines <- readLines(paste0(file, ".in"))
    for (i in seq_along(labels)) {
        mark <- which(lines %~% es("^(\\s*)/\\* MARK - #{labels[i]} \\*/"))
        if (length(mark) == 1) {
            padding <- groups()
            strings[[i]] <- paste0(ifelse(is.na(padding),"",padding), strings[[i]])
            lines <- c(lines[seq_len(mark-1)], strings[[i]], lines[-seq_len(mark)])
        }
    }
    writeLines(lines, file)
}

replacePlaceholder(file.path("src","main.cpp"), "Register callables", list(registerCallables))
replacePlaceholder(file.path("inst","include","RNiftiAPI.h"), c("Declare pointers", "Map pointers", "Define wrappers"), list(declarePointers, mapPointers, defineWrappers))
