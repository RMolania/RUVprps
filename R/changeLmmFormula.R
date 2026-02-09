#' Changes the formula for linear mixed model
#'
#' @param form Character. A character string specifying linear mixed model form.
#' @param out.put Character. A character string specifying the type of output of function. The options are `character`
#' and `sub.formula`. The default is set to `character`.
#' @param sub.set Character. A character string specifying which subset of the form should be extracted. The defulat is
#' set to `NULL`.

changeLmmFormula <- function(
        form,
        out.put = 'character',
        sub.set = NULL
        ){
    # Convert formula to character ####
    if (out.put == 'character'){
        formula.chara <- deparse(form)
        formula.chara <- gsub("[~()+|1]", "", formula.chara)
        formula.chara <- gsub("\\s+", " ", formula.chara)
        formula.chara <- strsplit(formula.chara, " ")[[1]]
        formula.chara <- formula.chara[formula.chara != ""]
        return(formula.chara)
    }
    if (out.put == 'sub.formula'){
        terms <- unlist(strsplit(paste(deparse(form), collapse = ""), "\\s*\\+\\s*"))
        terms <- gsub("~", "", terms)
        keep.terms <- terms[sapply(terms, function(term) any(sapply(sub.set, grepl, term)))]
        new.form <- paste("~", paste(keep.terms, collapse = " + "))
        return(new.form)
    }
}
