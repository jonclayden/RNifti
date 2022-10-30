#! /usr/bin/env Rscript

library(ore)

versions <- NULL

for (version in 1:2) {
    # Read header file containing function declarations
    lines <- readLines(file.path("inst", "include", "niftilib", es("nifti#{version}_io.h")))
    
    # Unwrap continuation lines
    wrapped <- which(lines %~% "^[^\\(]+\\)")
    while (length(wrapped) > 0) {
        lines[wrapped-1] <- paste(lines[wrapped-1], lines[wrapped])
        lines <- lines[-wrapped]
        wrapped <- which(lines %~% "^[^\\(]+\\)")
    }
    
    # Pattern-match each declaration, identifying groups corresponding to the
    # return type, function name and function arguments
    pieces <- groups(ore.search("^(\\w[\\w\\s]+?\\**)\\s*(\\w+)\\s*\\(([^\\)]+)\\)\\s*;\\s*(/\\*[^*]+\\*/\\s*)?$", lines))
    values <- ore.subst("^\\s+(.+?)\\s+", "\\1", pieces[,1,])
    values <- ore.subst("\\s+(\\s)", "\\1", values, all=TRUE)
    noreturn <- values == "void"
    symbols <- pieces[,2,]
    
    # Remove leading and trailing whitespace from the arg list, split at commas
    # and find the type and name of each argument
    args <- ore.subst("^\\s*(.+?)\\s*$", "\\1", pieces[,3,])
    args <- ore.subst("\\s+([\\s,])", "\\1", args, all=TRUE)
    argPieces <- groups(ore.search("(\\Avoid\\Z)|\\b([\\w\\s]+?\\**)\\s*(\\w+)\\s*(\\[\\d*\\])?(,|\\Z)", args, all=TRUE))
    
    # Create shortened names for the R-internal callable functions
    abbreviatedSymbols <- ore.subst("nifti", "nii", symbols)
    
    # Reconstruct the argument list: types only, names only and in full
    voidType <- sapply(argPieces, function(x) isTRUE(!is.na(x[,1])))
    argTypes <- ifelse(voidType, "void", NA_character_)
    argTypes[!voidType] <- sapply(argPieces[!voidType], function(x) paste0(x[,2], ifelse(!is.na(x[,4]),x[,4],""), collapse=", "))
    argNames <- sapply(argPieces, function(x) paste(x[,3],collapse=", "))
    args <- sapply(argPieces, function(x) paste0(x[,2], " ", x[,3], ifelse(!is.na(x[,4]),x[,4],""), collapse=", "))
    
    # Construct code fragments required for each function
    registerCallables <- paste0("R_RegisterCCallable(\"RNifti\", \"", abbreviatedSymbols, "\", (DL_FUNC) &", symbols, ");")
    declarePointers <- paste0("static ", values, "(*_", symbols, ")(", argTypes, ") = NULL;")
    mapPointers <- paste0("_", symbols, " = (", values, "(*)(", argTypes, ")) R_GetCCallable(\"RNifti\", \"", abbreviatedSymbols, "\");")
    
    # Construct wrapper functions
    defineWrappers <- rep(NA, length(argPieces))
    defineWrappers[voidType & !noreturn] <- paste0(values[voidType & !noreturn], " ", symbols[voidType & !noreturn], " (void) { NIFTILIB_WRAPPER_BODY_VOID(_", symbols[voidType & !noreturn], ") }")
    defineWrappers[!voidType & !noreturn] <- paste0(values[!voidType & !noreturn], " ", symbols[!voidType & !noreturn], " (", args[!voidType & !noreturn], ") { NIFTILIB_WRAPPER_BODY(_", symbols[!voidType & !noreturn], ", ", argNames[!voidType & !noreturn], ") }")
    defineWrappers[voidType & noreturn] <- paste0("void ", symbols[voidType & noreturn], " (void) { NIFTILIB_WRAPPER_BODY_VOID_NORETURN(_", symbols[voidType & noreturn], ") }")
    defineWrappers[!voidType & noreturn] <- paste0("void ", symbols[!voidType & noreturn], " (", args[!voidType & noreturn], ") { NIFTILIB_WRAPPER_BODY_NORETURN(_", symbols[!voidType & noreturn], ", ", argNames[!voidType & noreturn], ") }")
    defineWrappers <- na.omit(defineWrappers)
    
    versions <- c(versions, list(list(registerCallables=registerCallables, declarePointers=declarePointers, mapPointers=mapPointers, defineWrappers=defineWrappers)))
}

replacePlaceholder <- function (file, labels, vars, merge = FALSE) {
    lines <- readLines(paste0(file, ".in"))
    for (i in seq_along(labels)) {
        mark <- which(lines %~% es("^(\\s*)/\\* MARK - #{labels[i]} \\*/"))
        if (length(mark) == 1) {
            padding <- groups()
            strings <- lapply(1:2, function(version) paste0(ifelse(is.na(padding),"",padding), versions[[version]][[vars[i]]]))
            if (merge)
                lines <- c(lines[seq_len(mark-1)], unique(Reduce(c,strings)), lines[-seq_len(mark)])
            else
                lines <- c(lines[seq_len(mark-1)], "#if RNIFTI_NIFTILIB_VERSION == 1", strings[[1]], "#elif RNIFTI_NIFTILIB_VERSION == 2", strings[[2]], "#endif", lines[-seq_len(mark)])
        }
    }
    writeLines(lines, file)
}

replacePlaceholder(file.path("src","zzz.c"), "Register callables", "registerCallables", merge=TRUE)
replacePlaceholder(file.path("inst","include","RNiftiAPI.h"), c("Declare pointers", "Map pointers", "Define wrappers"), c("declarePointers", "mapPointers", "defineWrappers"))
