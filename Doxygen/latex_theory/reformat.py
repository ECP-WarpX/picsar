import re

# Open existing file
with open('theorycore.md') as f:
    text = f.read()

# Replace @ sign by \cite commands
text = text.replace( '@', '\cite ' )

# Replace symbol \bm which cannot be parsed by Doxygem
text = text.replace( r'\bm', r'\boldsymbol' )

# Replace inline math by proper math for Doxygen
text = text.replace( r'$$\begin{aligned}', r'\f{eqnarray*}{' )
text = text.replace( r'\end{aligned}$$', r'\f}' )
text = text.replace( '$', r'\f$' )

# Add header and footer for Doxygen
text = '/** @page "Theory"\n\n' + text + '\n*/\n'

# Write this to the page directory
with open('../pages/theory.txt', 'w') as f:
    f.write( text )
