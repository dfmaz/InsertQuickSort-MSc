# -*- coding: utf-8 -*-
"""
Práctica 2. Optimización y complejidad.
Daniel Fernández Martínez
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit

plt.style.use('default')
plt.rcParams.update({'font.size': 12})

" InsertionSort Algorithm"
def insertSort(array):
  
    # Traverse through 1 to len(arr)
    for i in range(1, len(array)):
  
        key = array[i]
  
        # Move elements of arr[0..i-1], that are
        # greater than key, to one position ahead
        # of their current position
        j = i-1
        while j >=0 and key < array[j] :
                array[j+1] = array[j]
                j -= 1
        array[j+1] = key
    return array

" QuickSort Algorithm"
def partition(arr, low, high):
    i = (low-1)         # index of smaller element
    pivot = arr[high]     # pivot
  
    for j in range(low, high):
  
        # If current element is smaller than or
        # equal to pivot
        if arr[j] <= pivot:
  
            # increment index of smaller element
            i = i+1
            arr[i], arr[j] = arr[j], arr[i]
  
    arr[i+1], arr[high] = arr[high], arr[i+1]
    return (i+1)
  
# The main function that implements QuickSort
# arr[] --> Array to be sorted,
# low  --> Starting index,
# high  --> Ending index
  
# Function to do Quick sort  
def quickSort(arr, low, high):
    if len(arr) == 1:
        return arr
    if low < high:
  
        # pi is partitioning index, arr[p] is now
        # at right place
        pi = partition(arr, low, high)
  
        # Separately sort elements before
        # partition and after partition
        quickSort(arr, low, pi-1)
        quickSort(arr, pi+1, high)

# Polynomial function for insertsort
def func1(N, a,b,c):
  return a + b*N + c*N**2

# Logarithmic function for quicksort
def func2(N, a,b,c):
  return a + b*N + c*N*np.log(N)

mode = 1

if mode == 0:

    # Study of the InsertionSort Algorithm
    """ Generate 1000 N-component vectors. Over these vectors, compute the CPU
     time needed by the previous two algorithms to order them. Compute the average
     of these times."""
    
    num_arr = 10000
    time_N = []
    k = 0
    for i in range(10, 150, 10):
        v = np.random.normal(0, 1, size=(num_arr, i)) # num_arr i-component vectors
    
        totalTime = []
        for j in range(len(v)):
            start = time.time()
            v_insertion = insertSort(v[j])
            end = time.time()
            totalTime.append(end-start)
        
        time_N.append(np.mean(totalTime))
        print('%s seconds needed for InsertionSort algorithm ' % time_N[k] + '(N = %d)'% i)
        k = k+1
    
    plt.plot(v_insertion)
    plt.title('Last vector sorted with InsertionSort')
    plt.xlabel('Vector index')
    plt.ylabel('Value')
    plt.tight_layout()
    plt.savefig('insertsorted.png', dpi=300)
    plt.show()
    
    """ Check that the average time obtained with the InsertSort algorithm 
    verifies the following law: t(N) = c1 + c2N + c3N^2 """
    x_data = np.linspace(10, 150, len(time_N))
    plt.plot(x_data, time_N, label='true curve')
    plt.title('InsertSort: $t(N)=c_1+c_2N+c_3N^2$')
    
    "Fit of time_N"
    # Model fit
    popt, pcov = curve_fit(func1, x_data, time_N)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    print(a, b, c)
    
    # Plot curve
    fit = func1(x_data,a,b,c)
    plt.plot(x_data, fit, label='fit')
    plt.xlabel('Components (N)')
    plt.ylabel('Time (s)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('insertfit.png', dpi=300)
    plt.show()
    
    "Errors"
    perr = np.sqrt(np.diag(pcov))/np.sqrt(len(x_data))
    nstd = 3
    popt_up = popt + nstd*perr
    popt_dw = popt - nstd*perr
    
    fit_up = func1(x_data, *popt_up)
    fit_dw = func1(x_data, *popt_dw)
    
    plt.plot(x_data, time_N, 'k', label='true curve')
    plt.plot(x_data, fit, 'r',  label='best fit')
    plt.fill_between(x_data, fit_up, fit_dw, alpha=.25, label='3-$\sigma$ interval')
    plt.title('InsertSort: $t(N)=c_1+c_2N+c_3N^2$')
    plt.xlabel('Components (N)')
    plt.ylabel('Time (s)')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig('insertfit_interval.png', dpi=300)
    plt.show()

if mode == 1:
    # Study of the QuickSort Algorithm
    # """ Generate 1000 N-component vectors. Over these vectors, compute the CPU
    #  time needed by the previous two algorithms to order them. Compute the average
    #  of these times."""
    num_arr = 10000
    time_N = []
    k = 0
    for i in range(10, 150, 10):
        v = np.random.normal(0, 1, size=(num_arr, i)) # num_arr i-component vectors
    
        totalTime = []
        for j in range(len(v)):
            start = time.time()
            quickSort(v[j], 0, i-1)
            end = time.time()
            totalTime.append(end-start)
        
        time_N.append(np.mean(totalTime))
        print('%s seconds needed for QuickSort algorithm ' % time_N[k] + '(N = %d)'% i)
        k = k+1
    
    plt.plot(v[i])
    plt.title('Last vector sorted with QuickSort')
    plt.xlabel('Vector index')
    plt.ylabel('Value')
    plt.tight_layout()
    plt.savefig('quicksorted.png', dpi=300)
    plt.show()
    
    """ Check that the average time obtained with the InsertSort algorithm 
    verifies the following law: t(N) = a + bN + cNlog(N)"""
    x_data = np.linspace(10, 150, len(time_N))
    plt.plot(x_data, time_N, label='true curve')
    plt.title('QuickSort: $t(N)=a + bN + cN \log{N}$')

    "Fit of time_N"
    # Model fit
    popt, pcov = curve_fit(func2, x_data, time_N)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    print(a, b, c)
    # Plot curve
    fit = func2(x_data,a,b,c)
    plt.plot(x_data, fit, label='fit')
    plt.xlabel('Components (N)')
    plt.ylabel('Time (s)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('quickfit.png', dpi=300)
    plt.show()
    
    "Errors"
    perr = np.sqrt(np.diag(pcov))/np.sqrt(len(x_data))
    nstd = 3
    popt_up = popt + nstd*perr
    popt_dw = popt - nstd*perr
    
    fit_up = func2(x_data, *popt_up)
    fit_dw = func2(x_data, *popt_dw)
    
    plt.plot(x_data, time_N, 'k', label='true curve')
    plt.plot(x_data, fit, 'r', label='best fit')
    plt.fill_between(x_data, fit_up, fit_dw, alpha=.25, label='3-$\sigma$ interval')
    plt.legend(loc='upper left')
    plt.title('QuickSort: $t(N)=a + bN + cN \log{N}$', )
    plt.xlabel('Components (N)')
    plt.ylabel('Time (s)')
    plt.tight_layout()
    plt.savefig('quickfit_interval.png', dpi=300)
    plt.show()