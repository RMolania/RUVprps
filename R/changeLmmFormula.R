#' Changes the formula for linear mixed model
#'
#' @param form TTT
#' @param out.put TTTT
#' @param sub.set TTTT

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
        terms <- unlist(strsplit(deparse(form), "\\s*\\+\\s*"))
        terms <- gsub("~", "", terms)
        keep.terms <- terms[sapply(terms, function(term) any(sapply(sub.set, grepl, term)))]
        new.form <- paste("~", paste(keep.terms, collapse = " + "))
        return(new.form)
    }
}

