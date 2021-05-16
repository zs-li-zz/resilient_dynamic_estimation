using PGFPlotsX
latexengine!(PGFPlotsX.PDFLATEX)
time_scale=101
i=3
time_axis=(0:time_scale-1)*0.02
figure = @pgf TikzPicture(
        Axis(
            PlotInc(
            { color => "black", },
            Table([ time_axis, Y[i,1:time_scale] ])
                ),
            PlotInc(
            { color => "black", },
            Table([ time_axis, Ya[i,1:time_scale] ])
                ),
            # PlotInc(
            # { color => "blue", },
            # Table([ time_axis, Xl[i,1:time_scale] ])
            #     ),
            PlotInc(
            { color => "blue", },
            Table([ time_axis, Ya[i,1:time_scale]-Y[i,1:time_scale] ])
                )
                )
                )

print_tex(figure)
# PGFPlotsX.pgfsave("cl_no.tex", figure)
