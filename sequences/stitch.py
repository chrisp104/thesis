import re

## FILES ##
counter = -1     # for identifying the CDR3s

abNames = open('abNames.txt', 'r')
heavyNames = open('heavyNames.txt', 'r')
lightNames = open('lightNames.txt', 'r')
ighv=open('IGHV.fasta','r')
iglkv=open('IGLKV.fasta','r')
cdrh3_file = open('CDRH3.txt', 'r')
cdrl3_file = open('CDRL3.txt', 'r')
output = open('output.txt', 'w')


## Loop through Ab names and stitch
for j in range(0, 75):

    # write the current Ab
    abName = abNames.next()
    output.write(abName)

    # Find heavy chain
    heavyName = heavyNames.next()
    heavyName = str(heavyName)
    heavyName = heavyName[0:len(heavyName)-1]
    ighv.seek(0)

    for heavy in ighv:
        if heavy.startswith('>'):
            info = heavy.split('|')
            name = info[1]
            if name == heavyName:
                chain = ighv.next()
                chain = chain[0:len(chain)-1]+ighv.next()

                index = 0

                # find where the CDR3 region starts in chain
                for i in range(1, len(chain)):
                    index = len(chain)-i
                    if (chain[index] == 'C'):
                        break

                output.write('H: '+chain[0:index])

                # append the CDRH3
                cdrh3 = cdrh3_file.next()
                output.write(cdrh3[0:len(cdrh3)-1])

                # append the CDRH3 generic tail
                output.write('GQGTLVTVSS\n')

                break


    # Find light chain
    lightName = lightNames.next()
    lightName = str(lightName)
    lightName = lightName[0:len(lightName)-1]
    iglkv.seek(0)

    for light in iglkv:
        if light.startswith('>'):
            info = light.split('|')
            name = info[1]
            if name == lightName:
                chain = iglkv.next()
                chain = chain[0:len(chain)-1]+iglkv.next()

                index = 0

                # find where the CDR3 region starts in chain
                for i in range(1, len(chain)):
                    index = len(chain)-i
                    if (chain[index] == 'C'):
                        break

                output.write('L: '+chain[0:index])

                # append the CDRL3
                cdrl3 = cdrl3_file.next()
                output.write(cdrl3[0:len(cdrl3)-1])

                # append the CDRL3 generic tail
                output.write('GQGTKVEIKRTV')

                break

    output.write('\n\n')

abNames.close()
heavyNames.close()
lightNames.close()
ighv.close()
iglkv.close()
cdrh3_file.close()
cdrl3_file.close()
output.close()

