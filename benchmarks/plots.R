library(dplyr)
library(ggplot2)
source('analysis.R')
source('multiplot.R')

res <- as_data_frame(results)
res_by_run <- group_by(res,progname,datafile,hashsize)
res_summary <- summarise(res_by_run, 
    cnstr_gmu=exp(mean(log(constructtime))),
    srch_gmu=exp(mean(log(searchtime))),
    mem_gmu=exp(mean(log(memory))),
    rt_gmu=exp(mean(log(runtime)))
)
for (i in c('DISTINCT.1mil','DISTINCT.500thou','SKEW.1mil')){
p1 <- ggplot(
	filter(res_summary,datafile == i), 
	aes(x=hashsize,y=rt_gmu, group=progname)) +
	geom_line(aes(colour=progname)) +
	geom_point(aes(colour=progname)) + 
	scale_x_log10(
	    breaks=2^(10:15),
	    labels=c(expression(2^10), expression(2^11),expression(2^12),
		expression(2^13), expression(2^14), expression(2^15) )
	) +
	labs(x='',y='',title='Run Time') +
	theme(
	    axis.text = element_text(colour = "black")
	)
p2 <- ggplot(
	filter(res_summary,datafile == i), 
	aes(x=hashsize,y=srch_gmu, group=progname)) +
	geom_line(aes(colour=progname)) +
	geom_point(aes(colour=progname)) + 
	scale_x_log10(
	    breaks=2^(10:15),
	    labels=c(expression(2^10), expression(2^11),expression(2^12),
		expression(2^13), expression(2^14), expression(2^15) )
	) +
	labs(x='',y='',title='Search Time') +
	theme(
	    axis.text = element_text(colour = "black")
	)
p3 <- ggplot(
	filter(res_summary,datafile == i), 
	aes(x=hashsize,y=cnstr_gmu, group=progname)) +
	geom_line(aes(colour=progname)) +
	geom_point(aes(colour=progname)) + 
	scale_x_log10(
	    breaks=2^(10:15),
	    labels=c(expression(2^10), expression(2^11),expression(2^12),
		expression(2^13), expression(2^14), expression(2^15) )
	) +
	labs(x='',y='',title='Construction Time') +
	theme(
	    axis.text = element_text(colour = "black")
	)

p4 <- ggplot(
	filter(res_summary,datafile == i), 
	aes(x=hashsize,y=mem_gmu, group=progname)) +
	geom_line(aes(colour=progname)) +
	geom_point(aes(colour=progname)) + 
	scale_x_log10(
	    breaks=2^(10:15),
	    labels=c(expression(2^10), expression(2^11),expression(2^12),
		expression(2^13), expression(2^14), expression(2^15) )
	) +
	labs(x='',y='Size in Mb',title='GC Memory') +
	theme(
	    axis.text = element_text(colour = "black")
	)

png(filename=paste(i,'.png',sep=''),width=600,height=500)
multiplot(p1,p2,p3,p4,cols=2)
dev.off()
}
