import numpy as np

# Create a small table that displays every bit value of a list of ifps and colors the background correspondingly
def ifpTable(ifps, header=False):
	cmap = ["FFFFFF", "FFA248", "007ADA", "A30092", "FF004F", "00C564", "5A7832", "00C0F7"]
	html = "<table border=1 width='50%'>"
	
	if header != False:
		html = html + "<tr>"
		for h in header:
			html = html  + "<td align='center'><b>" + str(h) + "</b></td>"
		html = html + "</tr>\n"

	for i in range(6,-1,-1):
		html = html + "<tr>"
		for ifp in ifps:
			bit = (ifp >> i) & 1
			html = html + "<td align='center' bgcolor='#" + cmap[(i+1)*bit] + "'>" + str(bit) + "</td>"
		html = html + "</tr>\n"

	html = html + "</table>\n"

	return html

def makeTableRow(content):
	html = "<tr>"
	for cell in content:
		cell = "<td align='center'>" + cell + "</td>"
		html = html + cell
	return html + "</tr>\n\r"

def includeImage(image, output):
	html = "<img src='" + image + "' style='width:50%;'></img>\n\r"
	output.write(html)

def generateHTMLOutput(output_path, html_matrix):
	#html_output = open("ifp.html", "w")
	html_output = open(output_path, "w")
	html_output.write("<html>\n\r<heading><title>IFP</title></heading>\n\r<font face='verdana'><body><center>\n\r")
	html_output.write("<h1>Interaction Fingerprint Viewer</h1>\n\r")
	includeImage("first_matrix_ever.png", html_output)
	html_output.write("<table border='0' width='80%' cellpadding='15'>\n\r")

	# Create a table 
	for i in range(1,100):
		html_output.write(makeTableRow(ifp_matrix, "FFFFFF"))

	html_output.write("</table>")
	html_output.write("</body></font>\n\r</html>")

	html_output.close()