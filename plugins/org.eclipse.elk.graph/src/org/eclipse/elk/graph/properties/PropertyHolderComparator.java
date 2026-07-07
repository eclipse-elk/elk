/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.properties;

import java.util.Comparator;

/**
 * A comparator for property holders that can uses a specific property for comparison.
 *
 * @author msp
 * @param <T> the type of the property used for comparison
 */
public final class PropertyHolderComparator<T extends Comparable<T>>
        implements Comparator<IPropertyHolder> {

    /**
     * Create a property holder comparator for the given property.
     * 
     * @param property a property
     * @param <U> the property type
     * @return a comparator that considers the given property
     */
    public static <U extends Comparable<U>> PropertyHolderComparator<U> with(
            final IProperty<U> property) {
        return new PropertyHolderComparator<U>(property);
    }
    
    /** the property used for comparing. */
    private final IProperty<T> property;
    
    /**
     * Create a property holder comparator for the given property.
     * 
     * @param property a property
     */
    private PropertyHolderComparator(final IProperty<T> property) {
        this.property = property;
    }
    
    @Override
    public int compare(final IPropertyHolder ph1, final IPropertyHolder ph2) {
        T p1 = ph1.getProperty(property);
        T p2 = ph2.getProperty(property);
        if (p1 != null && p2 != null) {
            return p1.compareTo(p2);
        } else if (p1 != null) {
            return -1;
        } else if (p2 != null) {
            return 1;
        } else {
            return 0;
        }
    }

}
