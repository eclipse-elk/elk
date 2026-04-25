/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.ArrayList;
import java.util.List;

/**
 * Base class for classes that want to randomly access a list of elements. Access to non-existent list elements is
 * explicitly permitted. The list is automatically enlarged by adding enough default values for the requested elements
 * to exist. Subclasses must implement {@link #provideDefault()} to provide these default values.
 * 
 * @param <T> type of elements in the list.
 */
public abstract class AbstractRandomListAccessor<T> {
    
    /** The list which contains all the elements. */
    private List<T> list;
    
    
    /**
     * Creates a new instance.
     */
    protected AbstractRandomListAccessor() {
        list = new ArrayList<>();
    }
    

    /**
     * Provides a default value for list items. This method is called whenever the list must be enlarged to be able
     * to access an element at a given index.
     * 
     * @return a default element.
     */
    protected abstract T provideDefault();
    
    /**
     * Returns the list item at the given index. This method is always guaranteed to return an element, provided that
     * the index is at least zero.
     * 
     * @param index the element's index.
     * @return the element at the given index.
     * @throws IndexOutOfBoundsException if the index is negative.
     */
    protected final T getListItem(final int index) {
        if (index < 0) {
            throw new IndexOutOfBoundsException("Invalid index: " + index);
        }
        
        ensureListSize(index + 1);
        return list.get(index);
    }
    
    /**
     * Sets the list element at the given index. If the index does not exist yet, the list is enlarged appropriately.
     * 
     * @param index the element's index.
     * @param value the new element.
     * @throws IndexOutOfBoundsException if the index is negative.
     */
    protected final void setListItem(final int index, final T value) {
        if (index < 0) {
            throw new IndexOutOfBoundsException("Invalid index: " + index);
        }
        
        if (index < list.size()) {
            list.set(index, value);
        } else {
            // Avoid creating a new default value which will be immediately overwritten anyway
            ensureListSize(index);
            list.add(value);
        }
    }
    
    /**
     * Returns the number of elements currently in the list.
     * 
     * @return number of elements.
     */
    protected final int getListSize() {
        return list.size();
    }
    
    /**
     * Removes all items from the list.
     */
    protected void clearList() {
        list.clear();
    }
    
    /**
     * Ensures that the list has the given size by adding enough default elements.
     * 
     * @param size the list's new size.
     */
    private void ensureListSize(final int size) {
        for (int i = list.size(); i < size; i++) {
            list.add(provideDefault());
        }
    }
    
}
