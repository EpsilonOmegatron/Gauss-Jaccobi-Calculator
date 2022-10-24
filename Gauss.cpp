#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main()
{   
    int size; //number of coefficients
    float *co_ef_swapper, sol_swapper, max_in_row, row_sum; //They carry row pointer, solution value respectively for swapping
    cout << "Enter number of co-efficients: ";
    cin >> size;
    
    //Build matrix of co-efficients
    float** coefficients = new float*[size];
    if(size)
    {
        for (int i = 0; i < size; i++)
        coefficients[i] = new float[size];
    }

    //Input co-efficients of one equation in one line seperated by spaces
    for (int i = 0; i < size; i++)
    {
        cout << "Enter co-efficients of equation " << i+1 << ": ";
        for (int j = 0; j < size; j++)
        {
            cin >> coefficients[i][j];
        }
        
    }

    //Build matrix of equation answers (insert them all in one line seperated by spaces)
    float solutions[size];
    cout << "Enter solutions in order: ";
    for(int i = 0; i < size; i++)
        {  
            cin >> solutions[i];
        }

    //Matrix validity check
    for(int i = 0; i < size; i++)
    {
        row_sum = 0;
        max_in_row = coefficients[i][0];
        
        for(int j = 0; j < size; j++)
        {
            row_sum = row_sum + abs(coefficients[i][j]);
            
            if(abs(coefficients[i][j]) > abs(max_in_row))
            {
                max_in_row = abs(coefficients[i][j]);
            }
        }

        if((row_sum - max_in_row) > max_in_row)
                {
                    cout << "\nThis system doesn't follow the jaccobi rule." << endl;
                    return 0;
                }
    }

    //Rearrange
    for (int j = 0; j < size; j++)
    {
        for (int i = 0; i < size; i++)
        {
            if(abs(coefficients[i][i]) < abs(coefficients[j][i]))
            {
                //Swap co-efficient row by pointers
                co_ef_swapper = coefficients[i];
                coefficients[i] = coefficients[j];
                coefficients[j] = co_ef_swapper;
                    
                //Swap solution by values
                sol_swapper = solutions[i];
                solutions[i] = solutions[j];
                solutions[j] = sol_swapper;
            }
        }
            
    }

    //Print matrix for value checking
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << setw(5) << coefficients[i][j] << "x" << j + 1;
        }
        cout << " = " << solutions[i] << endl;
    }

    //Gauss Jaccobi
    int iterations, iteration_counter = 0;
    float *iterationResult = new float [size];
    float *rowSum = new float [size];
    float *oldIteration = new float [size];
    
    cout << "Enter initial guess in order: ";
    for (int i = 0; i < size; i++)
    {
        cin >> iterationResult[i];
    }
    
    cout << "Enter number of iterations: ";
    cin >> iterations;
    
    while (iterations > 0)
    {
        cout << "Iteration " << ++iteration_counter << ":" << endl;
        
        for (int k = 0; k < size; k++)
        {
            oldIteration[k] = iterationResult[k]; //save past iteration in another array
            rowSum[k] = 0; //re-initialize sum for next iteration
        }

        for (int i = 0; i < size ; i++) //loop goes over rows
        {   
            for (int j = 0; j < size; j++) //loop goes over elements in a row
            {
                if(j != i) //for ignoring principal diagonal elements
                {rowSum[i] = rowSum[i] + (coefficients[i][j] * oldIteration[j]);} //multiply by iterationResult[j] for Seidel, oldIteration[j] for Jaccobi
            }
            
            iterationResult[i] = (1/coefficients[i][i]) * (solutions[i] - rowSum[i]); //result of iteration calculation with sum part included
            cout << "x" << i + 1 << " = " << iterationResult[i] << endl;
        }
        
        iterations--;
    }

    return 0;
}
