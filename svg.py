"""
用于从svg图片中提取路径信息
"""

from xml.dom import minidom

svg_file = "C:/Users/beta/Desktop/fourier/svg_files/浙大.svg"
doc = minidom.parse(svg_file)  # parseString also exists
path_strings = [path.getAttribute('d') for path
                in doc.getElementsByTagName('path')]
doc.unlink()
print(path_strings)
print(type(path_strings))
f=open("test.txt","w")
f.writelines(path_strings)
f.close()