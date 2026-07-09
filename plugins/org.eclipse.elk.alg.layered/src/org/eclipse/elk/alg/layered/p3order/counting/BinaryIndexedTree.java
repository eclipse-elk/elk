/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order.counting;

/**
 * Sorted list of integers storing values from 0 up to the maxNumber passed on creation. Adding, removing and indexOf
 * (and addAndIndexOf) is in O(log maxNumber).
 * <p/>
 * Implemented as a binary tree where each leaf stores the number of integers at the leaf index and each node stores the
 * number of values in the left branch of the node.
 *
 * @author alan
 */
public class BinaryIndexedTree {
    private int[] binarySums;
    private int[] numsPerIndex;
    private int size;
    private int maxNum;

    /**
     * Construct tree given maximum number of elements.
     *
     * @param maxNum
     *            maximum number elements.
     */
    public BinaryIndexedTree(final int maxNum) {
        this.maxNum = maxNum;
        binarySums = new int[maxNum + 1];
        numsPerIndex = new int[maxNum];
        size = 0;
    }

    /**
     * Increment given index.
     *
     * @param index
     *            The index to increment.
     */
    public void add(final int index) {
        size++;
        numsPerIndex[index]++;
        int i = index + 1;
        while (i < binarySums.length) {
            binarySums[i]++;
            i += i & -i;
        }
    }

    /**
     * Sum all entries before given index, i.e. index - 1.
     *
     * @param index
     *            Not included end index.
     * @return sum.
     */
    public int rank(final int index) {
        int i = index;
        int sum = 0;
        while (i > 0) {
            sum += binarySums[i];
            i -= i & -i;
        }
        return sum;
    }

    /**
     * @return size of tree.
     */
    public int size() {
        return size;
    }

    /**
     * Remove all entries for one index.
     *
     * @param index
     *            the index
     */
    public void removeAll(final int index) {
        int numEntries = numsPerIndex[index];
        if (numEntries == 0) {
            return;
        }
        numsPerIndex[index] = 0;
        size -= numEntries;
        int i = index + 1;
        while (i < binarySums.length) {
            binarySums[i] -= numEntries;
            i += i & -i;
        }
    }

    /**
     * Clears contents of tree.
     */
    public void clear() {
        binarySums = new int[maxNum + 1];
        numsPerIndex = new int[maxNum];
        size = 0;
    }

    /**
     * @return true if this tree is empty
     */
    public boolean isEmpty() {
        return size == 0;
    }
}
