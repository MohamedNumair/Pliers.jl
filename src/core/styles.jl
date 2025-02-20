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

function set_journal_theme(; fontsize = 12,      )

    journal_pub_theme = Theme(
    
    fontsize = fontsize,
    
    fonts = Attributes(
        :bold => Pliers.Makie.texfont(:bold),
        :bolditalic => Pliers.Makie.texfont(:bolditalic),
        :italic => Pliers.Makie.texfont(:italic),
        :regular => Pliers.Makie.texfont(:regular)
        ),

    
    Axis = (
        xticksmirrored = true,
        yticksmirrored = true,
        xtickalign=1,
        ytickalign=1,
        xminorgridvisible = true,
        yminorgridvisible = true,
    ),
    Scatter = (
        markersize = 10,

    ),

    Legend = (
        labelsize = 11,
        markersize = 4,
        margin = 1, 
        framevisible = true,
        colgap = 0,
        rowgap = -1,
        
        ),
)

set_theme!(journal_pub_theme)

end