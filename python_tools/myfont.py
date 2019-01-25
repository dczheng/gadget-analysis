#!/usr/bin/env python3


fonts = [
{'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 10,
},

]

label_font1 = fonts[0]


def main():
    import matplotlib
    matplotlib.use( 'agg' )
    import matplotlib.pyplot as plt
    import numpy as np

    print( fonts )

    x = np.linspace( -np.pi*2, np.pi*2, 100 )
    y = np.sin(x)
    plt.plot( x, y )

    plt.xlabel( 'x', fonts[0] )
    plt.ylabel( 'sin(x)', fonts[0] )
    plt.title( r'$sin(x)$', fonts[0] )

    plt.savefig( 'font_test.png' )

if __name__ == '__main__':
    main()
