### --- function: CQO Residual Plot 

cqo_resid_plot = function(x, smooth = T, points = T, legend = F ,logy = FALSE, logx = FALSE) {
        if (class(x) != "list") {
             x = list(x)   
        }
        
        require(data.table)
        require(VGAM)
        
        loop_list = list()
        for (i in seq_along(x)) {
                xx = x[[i]]
                xx_summary = summary(xx)
                env = as.matrix(xx_summary@x)
                C = concoef(xx)
                lv = numeric(length = nrow(env)) 
                for (j in 1:nrow(env)) lv[j] = t(C) %*% env[j,-1]      
                r = data.frame(resid(xx))   
                setDT(r)
                r[,lv:=lv]
                lv_id = which(names(r) == "lv")
                r = melt(r, id.vars = "lv", measure.vars = names(r)[-lv_id], variable.name = "taxon", value.name = "residual")
                r[,model_id := i]
                loop_list[[i]] = r
        }
        plot_data = rbindlist(loop_list)
        return_plot = 
                ggplot(data = plot_data, 
                       aes(x = lv, 
                           y = residual)
                       ) + 
                scale_color_viridis_d() + 
                xlab("linear predictor") + 
                geom_hline(yintercept = 0, linetype ="dashed") + 
                theme(legend.position = "none") + 
                facet_wrap(.~model_id)
        if (smooth) return_plot = return_plot + geom_smooth(aes(col = taxon), se = FALSE) 
        if (points) return_plot = return_plot + geom_point(aes(col = taxon), size = 2)
        if (legend) return_plot = return_plot + theme(legend.position = "right")
        if (logx)  return_plot = return_plot + scale_x_log10()
        return_plot
        
}
