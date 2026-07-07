/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later  
 *******************************************************************************/
package org.eclipse.elk.alg.topdownpacking;

import java.util.List;

/**
 * Provides an interface to store and access items in a grid structure. The grid size is bounded
 * by the size set. The index for the first column and row begins at 0. Empty cells of the grid
 * contain null.
 *
 */
public interface Grid<T> {
    
    /**
     * Stores an item at the requested position in the grid. If the requested index
     * is out of bounds an exception is thrown.
     * @param col index of requested column
     * @param row index of requested row
     * @param item object to be stored in the grid
     * @throws IndexOutOfBoundsException
     */
    public void put(int col, int row, T item) throws IndexOutOfBoundsException;
    
    /**
     * Returns the item stored at the requested position in the grid. If the requested index
     * is out of bounds an exception is thrown.
     * @param col index of requested column
     * @param row index of requested row
     * @return object stored at requested index
     * @throws IndexOutOfBoundsException
     */
    public T get(int col, int row) throws IndexOutOfBoundsException;
    
    /**
     * Returns the requested row of the grid as a {@link List}. If the requested row
     * is outside of the bounds of the grid an exception is thrown.
     * @param row the row number that should be returned
     * @return the list containing the entries of the row
     * @throws IndexOutOfBoundsException
     */
    public List<T> getRow(int row) throws IndexOutOfBoundsException;
    
    /**
     * Returns the requested column of the grid as a {@link List}. If the requested column
     * is outside of the bounds of the grid an exception is thrown.
     * @param col the column number that should be returned
     * @return the list containing the entries of the column
     * @throws IndexOutOfBoundsException
     */
    public List<T> getColumn(int col) throws IndexOutOfBoundsException;
    
    /**
     * Returns the number of columns in the grid.
     * @return the width of the grid
     */
    public int getColumns();
    
    /**
     * Returns the number of rows in the grid.
     * @return the height of the grid
     */
    public int getRows();
    
    /**
     * Defines the size of the grid. Access to coordinates outside the grid should throw an error.
     * 
     * @param cols the width of the grid
     * @param rows the height of the grid
     */
    public void setGridSize(int cols, int rows);

}
