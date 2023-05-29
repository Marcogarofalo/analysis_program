#define sorting_C


// A utility function to swap two elements
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}
/* This function takes last element as pivot, places
 *  the pivot element at its correct position in sorted
 *   array, and places all smaller (smaller than pivot)
 *  to left of pivot and all greater elements to right
 *  of pivot */
template<class T>
int partition (int *order, T*arr, int low, int high)
{
    double pivot = arr[order[high]];    // pivot
    int i = (low - 1);  // Index of smaller element
    
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[order[j]] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&order[i], &order[j]);
        }
    }
    swap(&order[i + 1], &order[high]);
    return (i + 1);
}

/* The main function that implements QuickSort
 * arr[] --> Array to be sorted,
 * low  --> Starting index,
 * high  --> Ending index */
template<class T>
void quickSort(int *order, T *arr, int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
         *          at right place */
        int pi = partition(order,arr, low, high);
        
        // Separately sort elements before
        // partition and after partition
        quickSort(order, arr, low, pi - 1);
        quickSort(order, arr, pi + 1, high);
    }
}

template void quickSort<int>(int *, int *, int, int);
template void quickSort<double>(int *, double *, int, int);