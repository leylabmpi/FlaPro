# Trimmed, project-local copy of leylabmpi/r_helper_lib :: functions_mb.R
# Vendored for FlaPro / notebooks/differential_analysis (IBD.ipynb).
# Upstream commit: 858d5b31f650cfea536205486ffb6cb1635a2f6e
# Contains ONLY the functions IBD.ipynb needs (4 of them); see UPSTREAM.md.
# Do not add functions here -- pull them from upstream and re-pin.

require(dplyr)
require(tibble)
require(tidyr)
require(purrr)
require(broom)
require(broom.mixed)   # tidy() for lmer fits
require(lme4)
require(lmerTest)
require(testthat)

# check the long table is full (includes the same numbers of zeros as if it were wide) before computing the prevalence
check_full = function(x, arg_feature_col) {
    tmp_n_samples = x %>% select(Sample) %>% distinct() %>% nrow()
    tmp_n_feats = x %>% select(arg_feature_col) %>% distinct() %>% nrow()
    test_that("Check the table is full before computing the prevalence", {    
        expect_equal(tmp_n_samples * tmp_n_feats,  x %>% nrow())
    })
}

adjust_via_residuals = function(arg_meta_samples, arg_features_long, arg_factor_to_adj_for, arg_feature_col = "feature", arg_value_col = "value") {
    check_full(arg_features_long, arg_feature_col)

    feat_df = 
        arg_features_long %>%          
        select(!!sym(arg_feature_col), Sample, !!sym(arg_value_col)) %>% 
        pivot_wider(names_from = !!sym(arg_feature_col), values_from = !!sym(arg_value_col), values_fill = 0) %>% 
        column_to_rownames("Sample") %>% 
        as.data.frame

    meta_df_tmp = data.frame(arg_meta_samples %>% column_to_rownames("Sample"))
    meta_df_tmp = meta_df_tmp[rownames(feat_df),]

    test_that("feat_df and arg_meta_samples have the same names", {
        expect_equal(rownames(feat_df), rownames(meta_df_tmp))
    })

    # collect the residuals
    lm_res_ = lm(reformulate(arg_factor_to_adj_for, response = "as.matrix(feat_df)"), meta_df_tmp)

    feat_df_adj = data.frame(resid(lm_res_))
    # non-alphanum characters appear to be replaced with dots by lm, so restore the original names
    colnames(feat_df_adj) = colnames(feat_df)

    tmp = feat_df_adj %>% rownames_to_column("Sample") %>%
        pivot_longer(cols = -Sample, names_to = arg_feature_col, values_to = arg_value_col) 

    tmp = tmp %>% inner_join(arg_features_long %>% select(-arg_value_col), by = c("Sample", arg_feature_col))   
}

do_lm_tidy = function(arg_lm_in, arg_formula, arg_response_col = "value", arg_feature_col = "feature") {    
    arg_lm_in %>%
        nest(data = -!!sym(arg_feature_col)) %>%    
        mutate(fit = map(data, function(x) {        
            lm1 = lm(reformulate(arg_formula, response = arg_response_col), data = x)
            return(lm1)
        }),
        tidied = map(fit, tidy)) %>%    
        select(!!sym(arg_feature_col), tidied) %>% 
        unnest(tidied) %>%
        filter(term != "(Intercept)") %>% 
        #filter(effect != "ran_pars") %>%            
        arrange(desc(term), p.value) #%>%
        #mutate(p.adj = p.adjust(p.value, method = "fdr"))   
}

do_lmer_tidy = function(arg_lm_in, arg_formula, arg_response_col = "value", arg_feature_col = "feature") {
    arg_lm_in %>%
        nest(data = -!!sym(arg_feature_col)) %>%    
        mutate(fit = map(data, function(x) {        
            lm1 = lmer(reformulate(arg_formula, response = arg_response_col), data = x)
            return(lm1)
        }),
        tidied = map(fit, tidy)) %>%    
        select(!!sym(arg_feature_col), tidied) %>% 
        unnest(tidied) %>%
        filter(term != "(Intercept)") %>% 
        filter(effect != "ran_pars") %>%            
        arrange(desc(term), p.value) #%>%
        #mutate(p.adj = p.adjust(p.value, method = "fdr"))   
}
