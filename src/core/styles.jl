# Printing Styles


function header(text::String)
    print(BOLD(UNDERLINE(BLUE_FG("$text\n"))))
end

function sub_header(text::String)
    print(BOLD(ITALICS(BLUE_FG("$text\n"))))
end

sub_sub_header(text::String) = print(BOLD(MAGENTA_FG("$text\n")))


function warning_text(message::String)
    print(ITALICS(YELLOW_FG("Warning: $message\n")))
end

function error_text(message::String)
    print(ITALICS(RED_FG("Error: $message\n")))
end
    

# PrettyTables Highlighting

highlight_row_label = Highlighter(
    f=(data, i, j) -> j == 1,
    crayon=Crayon(bold=true)
)

highlight_diagonal = Highlighter(
    f=(data, i, j) -> i == j-1,
    crayon=Crayon(foreground = :cyan, bold=true)
)

highlight_off_diagonal = Highlighter(
    f=(data, i, j) -> i != j-1,
    crayon=Crayon(foreground = :white, bold=true)
)

# for the results table in PMD
_highlight_results_status = Highlighter(
    f=(data, i, j) -> i==1 && string(data[i,j]) == "PF_CONVERGED",
    crayon=Crayon(foreground = :green, bold=true)
)


#=

 __       __            __        __           
/  \     /  |          /  |      /  |          
$$  \   /$$ |  ______  $$ |   __ $$/   ______  
$$$  \ /$$$ | /      \ $$ |  /  |/  | /      \ 
$$$$  /$$$$ | $$$$$$  |$$ |_/$$/ $$ |/$$$$$$  |
$$ $$ $$/$$ | /    $$ |$$   $$<  $$ |$$    $$ |
$$ |$$$/ $$ |/$$$$$$$ |$$$$$$  \ $$ |$$$$$$$$/ 
$$ | $/  $$ |$$    $$ |$$ | $$  |$$ |$$       |
$$/      $$/  $$$$$$$/ $$/   $$/ $$/  $$$$$$$/ 
                                               
                                               
                                               
=#


# Makie Themes

function set_journal_theme(; fontsize = nothing,      )

    inch = 96
    pt = 4/3
    cm = inch / 2.54
    ieeecolumn = 3.5 * inch
    ieee2column = 7.16 * inch

    # colors definition:

    # phase_red = Pliers.RGBAf(0.5,0.0,0.0,1.0)
    # phase_green = Pliers.RGBAf(0.0859375,0.3125,0.17578125,1.0)
    # phase_blue = Pliers.RGBAf(0.0,0.0,0.5,1.0)
    # neutral_black = Pliers.RGBAf(0.0,0.0,0.0,1.0)
    phase_red = :darkred
    phase_green = :darkgreen
    phase_blue = :darkblue
    neutral_black = :black
try
    global fontTermesRegular =  joinpath(dirname(dirname(pathof(Pliers))), "assets", "fonts", "gyre-opentype", "texgyretermes-regular.otf")
    global fontTermesBold =  joinpath(dirname(dirname(pathof(Pliers))), "assets", "fonts", "gyre-opentype", "texgyretermes-bold.otf")
    global fontTermesItalic =  joinpath(dirname(dirname(pathof(Pliers))), "assets", "fonts", "gyre-opentype", "texgyretermes-italic.otf")
    global fontTermesBoldItalic =  joinpath(dirname(dirname(pathof(Pliers))), "assets", "fonts", "gyre-opentype", "texgyretermes-bolditalic.otf")

catch
    global fontTermesRegular = Pliers.Makie.texfont(:regular)
    global fontTermesBold = Pliers.Makie.texfont(:bold)
    global fontTermesItalic = Pliers.Makie.texfont(:italic)
    global fontTermesBoldItalic = Pliers.Makie.texfont(:bolditalic)
    warning_text("Fonts not found for the IEEE fonts will default to ComputerModern")
end

    journal_pub_theme = Theme(
                                
                                fontsize = fontsize == nothing ? 8pt : fontsize,
                                palette = (color = [:darkred, :darkgreen, :darkblue, :black], marker = [:circle, :xcross]),
                                # fonts = Attributes(
                                #     :bold => Pliers.Makie.texfont(:bold),
                                #     :bolditalic => Pliers.Makie.texfont(:bolditalic),
                                #     :italic => Pliers.Makie.texfont(:italic),
                                #     :regular => Pliers.Makie.texfont(:regular)
                                #     ),
                                fonts = (; regular = fontTermesRegular, bold = fontTermesBold, italic = fontTermesItalic, bolditalic = fontTermesBoldItalic),

                                
                                Axis = (
                                    xticksmirrored = true,
                                    yticksmirrored = true,
                                    xtickalign=1,
                                    ytickalign=1,
                                    xminorgridvisible = true,
                                    yminorgridvisible = true,
                                ),
                                Scatter = (
                                    markersize = 5,
                                    cycle = [:color, :marker],
                                ),

                                Legend = (
                                    #labelsize = fontsize,
                                    markersize = 4,
                                    framevisible = true,
                                    colgap = -1,
                                    rowgap = -3,
                                    
                                    
                                    
                                    ),

                                
                                #palette = (color = [:red, :green, :blue, :black, :orange, :purple, :yellow, :cyan, :magenta, :gray]),
                                #Lines = (cycle = Cycle([:color, :linestyle], covary = true),)

                                #TODO: encode  the colors from Inkscape here 
)
theme_dark()
set_theme!(journal_pub_theme)


return inch, pt, cm, ieeecolumn, ieee2column

end


#=


IEEE Size
The size of a graphic refers to its dimensions (the width and height), which may be measured in inches, millimeters, or picas. Most charts, graphs, and tables are sized to be one column width or two columns width:

One column width: 3.5 inches, 88.9 millimeters, or 21 picas 
Two columns width: 7.16 inches, 182 millimeters, or 43 picas


=#

