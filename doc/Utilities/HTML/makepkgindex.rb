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
# Build an index.html file for all html files in the specified package.
#
# This script is used by doc/Developers/Modules/Makefile to build
# an index.html file for each package.
#

# Path to top of tree from Modules directory.
TREE_TOP = "../../../"

class Entry
  attr_reader :base, :file
  def initialize(file)
    @file = file
    @base = File.basename(file)
  end
  def <=>(e)
    @base <=> e.base
  end
  def get
    "<td><a href='#{@file}'>#{File.basename(@file, '.html')}</a></td>"
  end
end

class LatexEntry < Entry
  @@count = 0
  def initialize(file)
    super
    @@count += 1
  end
  def get
    "<td><a href='#{@file}'>#{File.basename(@file, '.html')}</a><sup>*</sup></td>"
  end
  def LatexEntry.count
    @@count
  end
end

def main

  packageName = ARGV[0]
  xmlDir = ARGV[1]
  texDir = ARGV[2]

  begin
    File.open("#{packageName}.html", "w") do |index|

      # Insert boilerplate at top of file
      index.print <<EndOfString
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
EndOfString

      version = File.open("../../edition.xml") do |f|
        /<edition>(.*)<\/edition>/.match(f.read())[1]
      end

      index.print("<br/>")

      # Generate index of modules.  Alphabetical order in 3 columns.
      entries = []
      Dir["#{xmlDir}/*.html"].each do |f|
          entries << Entry.new(f)
      end
		    
#       Dir["#{texDir}/*/*/*.html"].each do |f|
#         entries << LatexEntry.new(f)
#       end
      entries.sort!
      index.print("<table cellspacing='0' border='1'
cellpadding='2'> \n")
      # Generate page title and subtitle
      etable = []
      if entries.size() <= 25
	numCols = 1
	numRows = entries.size()
	numRowsArray = [ numRows ]
	etable[0] = entries
	extras = 0
      elsif entries.size() <= 50 or packageName == "Insight"
	numCols = 2
	numRows = entries.size() / numCols
	extras = entries.size() % numCols
	numRowsArray = [ numRows, numRows ]
	extras.times do |i|
	  numRowsArray[i] += 1
	end
	etable[0] = entries[0, numRowsArray[0]]
	etable[1] = entries[numRowsArray[0], numRowsArray[1]]
      else
	numCols = 3
	numRows = entries.size() / numCols
	extras = entries.size() % numCols
	numRowsArray = [ numRows, numRows, numRows ]
	extras.times do |i|
	  numRowsArray[i] += 1
	end
	etable = []
	etable[0] = entries[0, numRowsArray[0]]
	etable[1] = entries[numRowsArray[0], numRowsArray[1]]
	etable[2] = entries[numRowsArray[0] + numRowsArray[1], numRowsArray[2]]
      end

      index.print("
<thead>
<tr>
<th colspan=\"#{numCols}\">
<h1>#{packageName} Module Descriptions (v. #{version})</h1>
<h2><a class=\"alt\" href=\"index.html\">package index</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <a class=\"alt\" href=\"#{packageName}_bycat.html\">by category</a></h2>
</th>
</tr>
</thead>
<tbody>
\n")
      numRows.times do |i|
        index.print("<tr>\n")
        numCols.times do |j|
          index.print((etable[j][i]).get, "\n")
        end
        index.print("</tr>\n")
      end
      extras.times do |i|
	index.print("<tr>\n")
	index.print((etable[i][numRows]).get, "\n")
	index.print("</tr>\n")
      end

      index.print("</tbody>\n</table>\n")

      if LatexEntry.count > 0
        index.print("<p><sup>*</sup>Denotes additional documentation for a module.</p>\n")
      end

      # Generate end of file boilerplate.
      index.print <<EndOfString
<script type="text/javascript">doc.postContent();</script>
</body>
</html>
EndOfString

      end

  rescue => oops
    print("Unable to build #{packageName}.html for package\n")
    print("Reason: ", oops.message, "\n")
  end

end

main
