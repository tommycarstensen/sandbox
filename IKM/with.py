'''test2'''
x = 3

def main():

    with open('alice.txt','r') as fd:
        for line in fd:
            print line
    
    return

if __name__ == '__main__':
    main()
