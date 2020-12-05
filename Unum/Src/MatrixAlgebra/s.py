def run(filename):
    with open(filename, 'r') as file:
        biggest, l = 0, None
        for line in file:
            elements = line.split(',')
            if len(elements) > biggest:
                biggest=len(elements)
                l = elements
        print(l)

if __name__ == "__main__":
    run('A30_size1.csv')