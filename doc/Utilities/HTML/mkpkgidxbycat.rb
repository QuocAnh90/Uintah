#
#  The contents of this file are subject to the University of Utah Public
#  License (the "License"); you may not use this file except in compliance
#  with the License.
#
#  Software distributed under the License is distributed on an "AS IS"
#  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#  License for the specific language governing rights and limitations under
#  the License.
#
#  The Original Source Code is SCIRun, released March 12, 2001.
#
#  The Original Source Code was developed by the University of Utah.
#  Portions created by UNIVERSITY are Copyright (C) 2001, 1994
#  University of Utah. All Rights Reserved.
#

#
# Build an index of all html files in the specified package
# organized by category.
#
# This script is used by doc/Developers/Modules/Makefile to produce
# and index-by-cat.html file for each package.
#

# Path to top of tree from Modules directory.
TREE_TOP = "../../../"

# An Entry.  Consists of a value, is sortable, and provides the
# overridable entry function.
class Entry
  attr_reader :value
  def initialize(value)
    @value = value
  end
  def <=>(e)
    @value <=> e.value
  end
  def entry
    @value
  end
end

# An Entry for the name of a module category.
class CategoryEntry < Entry
  def initialize(category)
    super(category)
  end
  def entry
    "<td class=\"cattitle\">#{super()}</td>"
  end
end

# And Entry for a module description.
class ModuleEntry < Entry
  attr_reader :fileName
  def initialize(fileName)
    super(File.basename(fileName, ".html"))
    @fileName = fileName
  end
  def entry
    "<td><a href='#{@fileName}'>#{super()}</a></td>"
  end
  def category
    begin
      categoryP()
    rescue
      "REMOVE THESE!"
    end
  end
  def categoryP
    File.open(fileName.gsub(/\.html$/,".xml")) do |f|
      /category=\"(.*)\">/.match(f.read())[1]
    end
  end
end

# An Entry for a latex based module description.
class LatexEntry < ModuleEntry
  @@count = 0
  def initialize(fileName)
    super(fileName)
    @@count += 1
  end
  def entry
    "<td><a href='#{fileName()}'>#{File.basename(fileName, '.html')}</a><sup>*</sup></td>"
  end
  def categoryP
    xmlFile = (File.dirname(fileName) + "/../../../XML/" +
	       File.basename(fileName)).gsub(/\.html$/,".xml")
    File.open(xmlFile) do |f|
      /category=\"(.*)\">/.match(f.read())[1]
    end
  end
  def LatexEntry.count
    @@count
  end
end

def insertTopBoilerPlate(file, packageName)
  # Insert boilerplate at top of file
  file.print <<BOILER_PLATE_TEXT
<!--
The contents of this file are subject to the University of Utah Public
License (the "License"); you may not use this file except in compliance
with the License.

Software distributed under the License is distributed on an "AS IS"
basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
License for the specific language governing rights and limitations under
the License.

The Original Source Code is SCIRun, released March 12, 2001.

The Original Source Code was developed by the University of Utah.
Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
University of Utah. All Rights Reserved.
-->
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<script type="text/javascript" src="#{TREE_TOP}doc/Utilities/HTML/tools.js"></script>
<script type="text/javascript">var doc = new ModuleIndexDocument();</script>
<title>#{packageName} Module Descriptions</title>
</head>
<body>
<script type="text/javascript">doc.preContent();</script>
BOILER_PLATE_TEXT
end

def insertBottomBoilerPlate(file)
  # Generate end of file boilerplate.
  file.print <<BOILER_PLATE_TEXT
<script type="text/javascript">doc.postContent();</script>
</body>
</html>
BOILER_PLATE_TEXT
end

# Generate page of module description references in 1, 2, or 3
# columns.
def genPage(packageName, list)
  version = File.open("../../edition.xml") do |f|
    /<edition>(.*)<\/edition>/.match(f.read())[1]
  end
  File.open("#{packageName}_bycat.html", "w") do |file|
    insertTopBoilerPlate(file, packageName)

    file.print("<br/>\n")

    etable = []
    if list.size() <= 25
      numCols = 1
      numRows = list.size()
      numRowsArray = [ numRows ]
      etable[0] = list
      extras = 0
    elsif list.size() <= 50
      numCols = 2
      numRows = list.size() / numCols
      extras = list.size() % numCols
      numRowsArray = [ numRows, numRows ]
      extras.times do |i|
	numRowsArray[i] += 1
      end
      etable[0] = list[0, numRowsArray[0]]
      etable[1] = list[numRowsArray[0], numRowsArray[1]]
    else
      numCols = 3
      numRows = list.size() / numCols
      extras = list.size() % numCols
      numRowsArray = [ numRows, numRows, numRows ]
      extras.times do |i|
	numRowsArray[i] += 1
      end
      etable[0] = list[0, numRowsArray[0]]
      etable[1] = list[numRowsArray[0], numRowsArray[1]]
      etable[2] = list[numRowsArray[0] + numRowsArray[1], numRowsArray[2]]
    end

    file.print("<table cellspacing='0' border='1' cellpadding='2'>\n")
    # Generate page title and subtitle
    file.print("
<thead>
<tr>
<th colspan=\"#{numCols}\">
<h1>#{packageName} Module Descriptions (v. #{version})</h1>
<h2><a class=\"alt\" href=\"index.html\">package index</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <a class=\"alt\" href=\"#{packageName}.html\">alphabetically</a></h2>
</th>
</tr>
</thead>
<tbody>
")

    numRows.times do |r|
      file.print("<tr>\n")
      numCols.times do |c|
	file.print((etable[c][r]).entry, "\n")
      end
      file.print("</tr>\n")
    end
    extras.times do |i|
      file.print("<tr>\n")
      file.print((etable[i][numRows]).entry, "\n")
      file.print("</tr>\n")
    end

    file.print("</tbody>\n</table>\n")
    if LatexEntry.count > 0
      file.print("<p><sup>*</sup>Denotes additional documentation for a module.</p>\n")
    end

    insertBottomBoilerPlate(file)
  end
end

def main

  packageName = ARGV[0]
  xmlDir = ARGV[1]
  texDir = ARGV[2]

  begin
    # Collect xml descriptions into categories
    categories = {}
    Dir["#{xmlDir}/*.html"].each do |f|
      entry = ModuleEntry.new(f)
      if not categories.has_key?(entry.category)
	categories[entry.category] = []
      end
      categories[entry.category] << entry
    end
    
    # Add latex descriptions too.
    Dir["#{texDir}/*/*/*.html"].each do |f|
      entry = LatexEntry.new(f)
      if not categories.has_key?(entry.category)
	categories[entry.category] = []
      end
      categories[entry.category] << entry
    end

    # Sort modules in each category and insert category heading at the
    # front of each category list and create master list.
    master = []
    categories.each do |category,moduleList|
      moduleList.sort!
      moduleList.replace([ CategoryEntry.new(category) ] + moduleList)
    end

    # Now sort categories and make master module list.
    sortedCategories = categories.sort
    master = []
    sortedCategories.each do |category|
      master += category[1]
    end

    # Generate table of module descriptions.
    genPage(packageName, master)

  rescue => oops
    print("Unable to build #{packageName}.html for package\n")
    print("Reason: ", oops.message, "\n")
  end

end

main
