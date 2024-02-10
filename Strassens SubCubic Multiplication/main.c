//Strassens SubCubic Multiplication//
//Recursive new_size = size/2 for matrices of order higher than 2//
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MAX_SIZE 10

void matrix_multiplication(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size);
void matrix_addition(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size);
void matrix_subtract(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size);
void strassen(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size);
void display_matrix(int matrix[MAX_SIZE][MAX_SIZE],int size);

//ADDITION OF MATRICES//
void matrix_addition(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            c[i][j]=a[i][j]+b[i][j];
        }
    }
}

//Subtraction of Matrices//
void matrix_subtract(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            c[i][j]=a[i][j]-b[i][j];
        }
    }
}

void matrix_multiplication(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size)
{  if(size<=2)
 {
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            c[i][j]=0;

            for(int k=0;k<size;k++)
            {
                c[i][j]+=a[i][k]*b[k][j];
            }
        }
    }
 }
 else
    {
     strassen(a,b,c,size);
    }
}

void strassen(int a[MAX_SIZE][MAX_SIZE], int b[MAX_SIZE][MAX_SIZE], int c[MAX_SIZE][MAX_SIZE],int size)
{
    if(size<=2)
    {
        matrix_multiplication(a,b,c,size);
        return;
    }

    int new_size=size/2;

    //Submatrices of A//
    int a11[MAX_SIZE][MAX_SIZE], a12[MAX_SIZE][MAX_SIZE], a21[MAX_SIZE][MAX_SIZE], a22[MAX_SIZE][MAX_SIZE];

    //SubMatrices of B//
    int b11[MAX_SIZE][MAX_SIZE], b12[MAX_SIZE][MAX_SIZE], b21[MAX_SIZE][MAX_SIZE], b22[MAX_SIZE][MAX_SIZE];


    //Submatrices of C//
    int c11[MAX_SIZE][MAX_SIZE], c12[MAX_SIZE][MAX_SIZE], c21[MAX_SIZE][MAX_SIZE], c22[MAX_SIZE][MAX_SIZE];



    //Temporary Matrices for intermediate results//
    int temp1[MAX_SIZE], temp2[MAX_SIZE][MAX_SIZE], temp3[MAX_SIZE][MAX_SIZE], tempP5[MAX_SIZE][MAX_SIZE], tempP6[MAX_SIZE][MAX_SIZE], tempP7[MAX_SIZE][MAX_SIZE];

    //Partioning the matrices a and b//
    for(int i=0;i<new_size;i++)
    {
        for(int j=0;j<new_size;j++)
        {
            a11[i][j]=a[i][j];
            a12[i][j]=a[i][j+new_size];
            a21[i][j]=a[i+new_size][j];
            a22[i][j]=a[i+new_size][j+new_size];

            b11[i][j]=b[i][j];
            b12[i][j]=b[i][j+new_size];
            b21[i][j]=b[i+new_size][j];
            b22[i][j]=b[i+new_size][j+new_size];
        }
    }


    //P1:(A11+A22)*(B11+B22)//
    matrix_addition(a11,a22,temp1,new_size);
    matrix_addition(b11,b22,temp2,new_size);
    strassen(temp1,temp2,c11,new_size);

    //P2:(A21+A22)*(B11)
    matrix_addition(a21,a22,temp1,new_size);
    strassen(temp1,b11,c12,new_size);

    //P3:(A11)*(B12-B22)
    matrix_subtract(b12,b22,temp1,new_size);
    strassen(a11,temp1,c21,new_size);

    //P4:(A22)*(B21-B11)
    matrix_subtract(b21,b11,temp1,new_size);
    strassen(a22,temp1,c22,new_size);

    //P5:(A11+A12)*(B22)
    matrix_addition(a11,a12,temp1,new_size);
    strassen(temp1,b22,tempP5,new_size);

    //P6:(A21-A11)*(B11+B12)
    matrix_subtract(a21,a11,temp1,new_size);
    matrix_addition(b11,b12,temp2,new_size);
    strassen(temp1,temp2,tempP6,new_size);

    //P7:(A12-A22)*(B21+B22)
    matrix_subtract(a12,a22,temp1,new_size);
    matrix_addition(b21,b22,temp2,new_size);
    strassen(temp1,temp2,tempP7,new_size);

    //C11:P1+P4+P7-P5
    matrix_addition(c11,c22,temp1,new_size);
    matrix_addition(temp1,tempP7,temp2,new_size);
    matrix_subtract(temp2,tempP5,c11,new_size);

    //C12:P3+P5
    matrix_addition(c21,tempP5,c12,new_size);

    //C21:P2+P4
    matrix_addition(c12,c22,c21,new_size);

    //C22:P1+P3+P6-P2 */
    matrix_addition(c11,c21,temp1,new_size);
    matrix_addition(temp1,tempP6,temp2,new_size);
    matrix_subtract(temp2,c12,c22,new_size);

    //Recombining the C Matrix//
    for(int i=0;i<new_size;i++)
    {
        for(int j=0;j<new_size;j++)
    {
    c[i][j]=c11[i][j];
    c[i][j+new_size]=c12[i][j];
    c[i+new_size][j]=c21[i][j];
    c[i+new_size][j+new_size]=c22[i][j];
    }
    }

}


void display_matrix(int matrix[MAX_SIZE][MAX_SIZE], int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%d\t", matrix[i][j]); // Print each element followed by a tab character
        }
        printf("\n"); // Move to the next row after printing each row
    }
}






int main()
{
    int a[MAX_SIZE][MAX_SIZE],b[MAX_SIZE][MAX_SIZE],c[MAX_SIZE][MAX_SIZE];
    int size;

    printf("Enter the size of the matrices:");
    scanf("%d",&size);

    printf("Enter the Elements of the First Matrix:");
    printf("\n");
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            scanf("%d",&a[i][j]);
        }
    }

    printf("Enter the Elements of the Second Matrix:");
    printf("\n");

    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            scanf("%d",&b[i][j]);
        }
    }

    matrix_multiplication(a,b,c,size);

    printf("Resultant Matrix After Multiplication:");
    printf("\n");
    display_matrix(c,size);

    return 0;
}




