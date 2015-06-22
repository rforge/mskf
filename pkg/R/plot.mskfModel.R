plot.mskfModel <- function(x, cex = 3.5, length=0.1, col = 1:mdl$nm) {
	mdl = x
	plot.new()
	plot.window(c(-2,1), c(-1.25,1))


	if (mdl$ny > 0) {
		y.locx = seq(0, 1, by = 1 / (mdl$ny + 1));
		y.locx = y.locx[-c(1,length(y.locx))]
		y.locy = rep(0, mdl$ny);
		points(y.locx, y.locy, pch = 0, cex=cex)
		points(-y.locx, y.locy, pch = 0, cex=cex)
	}

	if (mdl$ne > 0) {
		eta.locx = rep(0.5, mdl$ne)
		eta.locy = seq(0, 1, by = 1 / (mdl$ne + 1)) + 0.5 / (mdl$ne + 1)
		eta.locy = eta.locy[-c(1,length(eta.locy))]
		points(eta.locx, eta.locy, pch=1, cex=cex)
		points(-eta.locx, eta.locy, pch=1, cex=cex)
	}

	if (mdl$nx > 0) {
		x.locx = rep(0.5, mdl$nx)
		x.locy = seq(-1, 0, by = 1 / (mdl$nx + 1)) - 0.5 / (mdl$ne + 1)
		x.locy = x.locy[-c(1,length(x.locy))]
		points(x.locx, x.locy, pch=5, cex=cex)
		points(-x.locx, x.locy, pch=5, cex=cex)
	}

	dr = seq(-1, 1, len=mdl$nm) * ifelse(length(diff(y.locx)) > 0, diff(y.locx)[1]/15, 0.0125)  # regime offsets
	de = seq(-1, 1, len=mdl$ne) * ifelse(length(diff(y.locx)) > 0, diff(y.locx)[1]/13,  0) # eta offsets

	for (fun in c(maW, paW)) {
		W = maW(mdl);
		lty = identical(fun, paW) + 1;
		if (!is.null(W)) {
			for (i in 1:mdl$ny)
				for (j in 1:mdl$ne)
					for (k in 1:mdl$nm)
						if(W[i,j,k] != 0) {
							o = c(x = eta.locx[j]+de[j]+dr[k], y = eta.locy[j]);
							e = c(x = y.locx[i]+de[j]+dr[k], y = y.locy[i]);
							x = 'x'; y = 'y';
							o = o + 0.15 * (e-o);
							arrows(o[x], o[y], e[x], e[y]+.1, col=col[k], lty=lty, length=length);
							o = c(x = sort(-eta.locx)[j]+de[j]+dr[k], y = eta.locy[j]);
							e = c(x = sort(-y.locx)[i]+de[j]+dr[k], y = y.locy[i]);
							o = o + 0.15 * (e-o);
							arrows(o[x], o[y], e[x], e[y]+.1, col=col[k], lty=lty, length=length);
						}
		}
	}


	for (fun in c(maB, paB)) {
		B = fun(mdl);
		lty = identical(fun, paB) + 1;
		if (!is.null(B)) {
			for (i in 1:mdl$ny)
				for (j in 1:mdl$nx)
					for (k in 1:mdl$nm)
						if(B[i,j,k] != 0) {
							o = c(x = x.locx[j]+de[j]+dr[k], y = x.locy[j]);
							e = c(x = y.locx[i]+de[j]+dr[k], y = y.locy[i]);
							x = 'x'; y = 'y';
							o = o + 0.1 * (e-o);
							arrows(o[x], o[y], e[x], e[y]-.1, col=col[k], lty=lty, length=length);
							o = c(x = sort(-x.locx)[j]+de[j]+dr[k], y = x.locy[j]);
							e = c(x = sort(-y.locx)[i]+de[j]+dr[k], y = y.locy[i]);
							o = o + 0.1 * (e-o)
							arrows(o[x], o[y], e[x], e[y]-.1, col=col[k], lty=lty, length=length);
						}
		}
	}

	for (fun in c(maH, paH)) {
		H = fun(mdl);
		lty = identical(fun, paH) + 1
		if (!is.null(H)) {
			for (i in 1:mdl$ne)
				for (j in 1:mdl$ne)
					for (k in 1:mdl$nm)
						if(H[i,j,k] != 0) {
							o = c(x = eta.locx[j], y = eta.locy[j]+de[i]+dr[k]);
							e = c(x = sort(-eta.locx)[i], y = eta.locy[i]+de[i]+dr[k]);
							x = 'x'; y = 'y';
							o = o + 0.1 * (e-o);
							arrows(o[x], o[y], e[x]+.1, e[y], col=col[k], lty=lty, length=length);
						}
		}
	}

	text(-1.5, eta.locy, parse(text = paste("eta[",1:mdl$ne,"]",sep="")))
	text(-1.5,   y.locy, parse(text = paste(paste("y[",if(mdl$ny>2) c(1,mdl$ny) else 1:mdl$ny,"]",sep=""),collapse=if(mdl$ny>2) "~...~" else "~")))
	if (mdl$nx > 0) {
		text(-1.5,   x.locy, parse(text = paste("x[",1:mdl$nx,"]",sep="")))
		text(-eta.locx, -1.25, expression(t))
		text( eta.locx, -1.25, expression(t-1))
	}
	else {
		text(-eta.locx, -0.5, expression(t))
		text( eta.locx, -0.5, expression(t-1))
	}

	legend(-2, -0.9, paste("regime", 1:mdl$nm), fill = 1:mdl$nm, cex=.7);

}
