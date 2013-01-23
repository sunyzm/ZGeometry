import sys

def main():
    if len(sys.argv) != 2:
        return
    else: 
        filename = sys.argv[1]
        file = open(filename, 'r')
        strold = file.read()
        file.close()
        open(filename,'w').write(strold.replace(' ','\n'))

if __name__ == '__main__':
    main()
