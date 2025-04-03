using PlutoSliderServer

# Export notebooks to HTML
practice_notebooks_dir = "practices"
export_dir = "public"

# Create export directory and .nojekyll file
mkpath(export_dir)
open(joinpath(export_dir, ".nojekyll"), "w") do io end

# Export all .jl notebooks
PlutoSliderServer.export_directory(practice_notebooks_dir, Export_output_dir = export_dir)

# Generate index.html with links
index_content = """
<!DOCTYPE html>
<html>
<head>
    <title>Задачи и материалы по курсу "Методы исследования динамических систем"</title>
    <style>
        body { font-family: sans-serif; margin: 2em; }
        a { color: #0366d6; }
    </style>
</head>
<body>
    <h1>Практические задания</h1>
    <ul>
"""

for file in readdir(export_dir)
    if endswith(file, ".html") && file ≠ "index.html"
		global index_content
	    index_content *= """<li><a href="$file">$(replace(file, ".html" => ""))</a></li>\n"""
    end
end

index_content *= """
    </ul>
</body>
</html>
"""

write(joinpath(export_dir, "index.html"), index_content)
