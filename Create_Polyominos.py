import os
import sys
import math
import itertools
import numpy as np
import matplotlib.pyplot as plt

num_of_squares_list = np.linspace(1, 15, 15, dtype=int)

for num_of_squares in num_of_squares_list:
        num_of_polyominos = 0

        max_rows = math.ceil(num_of_squares / 2)

        matrix_sizes = []
        polyominos = []

        solution_folder = str(num_of_squares) + '_square_pieces'
        if not os.path.exists(solution_folder):
                os.makedirs(solution_folder)

        pieces_per_batch = 1000

        sub_solution_folder_png = solution_folder + '/Pictures_1_to_' + str(pieces_per_batch)
        if not os.path.exists(sub_solution_folder_png):
                os.makedirs(sub_solution_folder_png)

        for i in range(1, max_rows + 1, 1):
                for j in range(1, num_of_squares + 1 - i + 1, 1):
                        if i * j < num_of_squares or i * j > max_rows * (num_of_squares - max_rows + 1):
                                continue

                        if [j, i] not in matrix_sizes:
                                matrix_sizes.append([i, j])
                        else:
                                continue

                        print("Checking [{}, {}] shapes...".format(i, j))

                        sub_solution_folder_ini = solution_folder + '/Pieces'
                        if not os.path.exists(sub_solution_folder_ini):
                                os.makedirs(sub_solution_folder_ini)
                        
                        file = open(sub_solution_folder_ini + '/Pieces_' + str(i) + 'x' + str(j) + '_1.ini', 'w')
                        counter = 1

                        cell_numbers = np.linspace(1, i * j, i * j, dtype=int)
                        polyominos = []
                        for combination in itertools.combinations(cell_numbers, num_of_squares):
                                polyominos.append(np.zeros([i, j], dtype=int))

                                # Create the pieces
                                for c in combination:
                                        row = int((c-1) / j)
                                        col = (c-1) % j

                                        polyominos[-1][row][col] = 1
                
                                # Remove the pieces if it has column of zeros
                                if (~polyominos[-1].any(axis=0)).any():
                                        del polyominos[-1]
                                        continue
                                
                                # Remove the pieces if it has row of zeros
                                if (~polyominos[-1].any(axis=1)).any():
                                        del polyominos[-1]
                                        continue

                                # Remove disjoint pieces
                                r = 0
                                c = 0
                                while polyominos[-1][r][c] != 1:
                                        c = c + 1
                                        if c == len(polyominos[-1][0]):
                                                c = 0
                                                r = r + 1

                                total_ones = 1
                                queue = []
                                checked_queue = []
                                queue.append([r, c])
                                checked_queue.append([r, c])
                                while queue != []:
                                        if queue[0][0] + 1 < len(polyominos[-1]) and polyominos[-1][queue[0][0] + 1][queue[0][1]] == 1 and [queue[0][0] + 1, queue[0][1]] not in checked_queue:
                                                total_ones = total_ones + 1
                                                queue.append([queue[0][0] + 1, queue[0][1]])
                                                checked_queue.append([queue[0][0] + 1, queue[0][1]])
                                        
                                        if queue[0][1] + 1 < len(polyominos[-1][0]) and polyominos[-1][queue[0][0]][queue[0][1] + 1] == 1 and [queue[0][0], queue[0][1] + 1] not in checked_queue:
                                                total_ones = total_ones + 1
                                                queue.append([queue[0][0], queue[0][1] + 1])
                                                checked_queue.append([queue[0][0], queue[0][1] + 1])
                                        
                                        if queue[0][0] - 1 > -1 and polyominos[-1][queue[0][0] - 1][queue[0][1]] == 1 and [queue[0][0] - 1, queue[0][1]] not in checked_queue:
                                                total_ones = total_ones + 1
                                                queue.append([queue[0][0] - 1, queue[0][1]])
                                                checked_queue.append([queue[0][0] - 1, queue[0][1]])
                                        
                                        if queue[0][1] - 1 > -1 and polyominos[-1][queue[0][0]][queue[0][1] - 1] == 1 and [queue[0][0], queue[0][1] - 1] not in checked_queue:
                                                total_ones = total_ones + 1
                                                queue.append([queue[0][0], queue[0][1] - 1])
                                                checked_queue.append([queue[0][0], queue[0][1] - 1])

                                        queue.pop(0)

                                if total_ones < num_of_squares:
                                        del polyominos[-1]
                                        continue

                                # Remove duplicates
                                continuer = False
                                # Spin the piece 3 times
                                check_polyominos = np.zeros([i, j])
                                check_polyominos = polyominos[-1]
                                for s in range(3):
                                        check_polyominos = np.transpose(check_polyominos)
                                        check_polyominos = check_polyominos[::-1]
                                        if any(np.array_equal(check_polyominos, i) for i in polyominos[:-1]):
                                                del polyominos[-1]
                                                continuer = True
                                                break
                                if continuer:
                                        continue

                                # Flip the piece and spin it again 3 times
                                check_polyominos = polyominos[-1][::-1]
                                if any(np.array_equal(check_polyominos, i) for i in polyominos[:-1]):
                                        del polyominos[-1]
                                        continue
                                for s in range(3):
                                        check_polyominos = np.transpose(check_polyominos)
                                        check_polyominos = check_polyominos[::-1]
                                        if any(np.array_equal(check_polyominos, i) for i in polyominos[:-1]):
                                                del polyominos[-1]
                                                continuer = True
                                                break
                                if continuer:
                                        continue

                                if len(polyominos) > 1:
                                        for s in range(i+1):
                                                sys.stdout.write("\033[F")
                                print("New shape found:")
                                print(polyominos[-1])

                                # Save the pieces (ini)
                                file.write('[Piece_' + str(counter) + ']\n')
                                for s in range(len(polyominos[-1])):
                                        file.write('row' + str(s) + ' = ')
                                        for ss in range(len(polyominos[-1][s])):
                                                file.write(str(polyominos[-1][s][ss]) + ' ')
                                        file.write('\n')
                                file.write('\n')

                                if counter % pieces_per_batch == 0:
                                        file.close()
                                        file = open(sub_solution_folder_ini + '/Pieces_' + str(i) + 'x' + str(j) + '_' + str(int(counter / pieces_per_batch) + 1) + '.ini', 'w')

                                # Save the pieces (png)
                                num_of_polyominos = num_of_polyominos + 1

                                if (num_of_polyominos > 1) and (num_of_polyominos - 1) % pieces_per_batch == 0:
                                        sub_solution_folder_png = solution_folder + '/Pictures_' + str(num_of_polyominos) + '_to_' + str(num_of_polyominos + pieces_per_batch - 1)
                                        if not os.path.exists(sub_solution_folder_png):
                                                os.makedirs(sub_solution_folder_png)

                                plt.figure(figsize=(j, i))
                                plt.xlim(0, j)
                                plt.ylim(0, i)
                                plt.gca().set_aspect('equal')
                                plt.axis('off')
                                
                                for ii in range(i):
                                        for jj in range(j):
                                                if ii == 0:
                                                        if polyominos[-1][ii][jj] == 1:
                                                                xs = [jj, jj + 1]
                                                                ys = [i - ii, i - ii]
                                                                plt.plot(xs, ys, color = 'black')
                                                
                                                if ii == i - 1:
                                                        if polyominos[-1][ii][jj] == 1:
                                                                xs = [jj, jj + 1]
                                                                ys = [i - ii - 1, i - ii - 1]
                                                                plt.plot(xs, ys, color = 'black')
                                                
                                                if ii + 1 < i:
                                                        if (polyominos[-1][ii][jj] == 1 and polyominos[-1][ii+1][jj] == 1) or \
                                                        (polyominos[-1][ii][jj] == 1 and polyominos[-1][ii+1][jj] == 0) or \
                                                        (polyominos[-1][ii][jj] == 0 and polyominos[-1][ii+1][jj] == 1):
                                                                xs = [jj, jj + 1]
                                                                ys = [i - ii - 1, i - ii - 1]
                                                                plt.plot(xs, ys, color = 'black')
                                                
                                                if jj == 0:
                                                        if polyominos[-1][ii][jj] == 1:
                                                                xs = [jj, jj]
                                                                ys = [i - ii, i - ii - 1]
                                                                plt.plot(xs, ys, color = 'black')
                                                
                                                if jj == j - 1:
                                                        if polyominos[-1][ii][jj] == 1:
                                                                xs = [jj + 1, jj + 1]
                                                                ys = [i - ii, i - ii - 1]
                                                                plt.plot(xs, ys, color = 'black')
                                                
                                                if jj + 1 < j:
                                                        if (polyominos[-1][ii][jj] == 1 and polyominos[-1][ii][jj+1] == 1) or \
                                                        (polyominos[-1][ii][jj] == 1 and polyominos[-1][ii][jj+1] == 0) or \
                                                        (polyominos[-1][ii][jj] == 0 and polyominos[-1][ii][jj+1] == 1):
                                                                xs = [jj + 1, jj + 1]
                                                                ys = [i - ii, i - ii - 1]
                                                                plt.plot(xs, ys, color = 'black')
                                                
                                                if polyominos[-1][ii][jj] == 1:
                                                        rectangle = plt.Rectangle((jj, i - ii - 1), 1, 1, fc='green')
                                                        plt.gca().add_patch(rectangle)

                                plt.savefig(sub_solution_folder_png + '/Piece_' + str(num_of_polyominos) + '.png')
                                plt.close()

                                counter = counter + 1

                        file.close()

                        print("Found {} new shapes of size [{}, {}], and total number of {} shapes.".format(len(polyominos), i, j, num_of_polyominos))

        print("Total number of polyominos that found with {} sqaures is: {}".format(num_of_squares, num_of_polyominos))

