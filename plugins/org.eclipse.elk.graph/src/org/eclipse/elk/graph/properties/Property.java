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

import org.eclipse.elk.graph.util.ElkReflect;

/**
 * Definition of a property that can be set on an {@link IPropertyHolder}. A property is described
 * by its identifier, which is a String. Properties may have a default value which is returned by
 * a property holder if the property is not set on it. Properties may also define lower and upper
 * bounds that restrict which values are considered legal.
 *
 * @param <T> type of the property
 * @author msp
 */
public class Property<T> implements IProperty<T>, Comparable<IProperty<?>> {
    
    /** identifier of this property. */
    private final String id;
    /** the default value of this property. */
    private T defaultValue;
    /** the lower bound of this property. */
    private Comparable<? super T> lowerBound = NEGATIVE_INFINITY;
    /** the upper bound of this property. */
    private Comparable<? super T> upperBound = POSITIVE_INFINITY;
    
    /**
     * Creates a property with given identifier and {@code null} as default value.
     * 
     * @param theid the identifier
     */
    public Property(final String theid) {
        this.id = theid;
    }
    
    /**
     * Creates a property with given identifier and default value.
     * 
     * @param theid the identifier
     * @param thedefaultValue the default value
     */
    public Property(final String theid, final T thedefaultValue) {
        this(theid);
        this.defaultValue = thedefaultValue;
    }
    
    /**
     * Creates a property using another property as identifier, but changing
     * the default value.
     * 
     * @param other another property
     * @param thedefaultValue the new default value
     */
    public Property(final IProperty<T> other, final T thedefaultValue) {
        this(other.getId(), thedefaultValue);
    }
    
    /**
     * Creates a property with given identifier, default value, and lower bound.
     * 
     * @param theid the identifier
     * @param thedefaultValue the default value
     * @param thelowerBound the lower bound
     */
    public Property(final String theid, final T thedefaultValue,
            final Comparable<? super T> thelowerBound) {
        this(theid, thedefaultValue);
        if (thelowerBound != null) {
            this.lowerBound = thelowerBound;
        }
    }
    
    /**
     * Creates a property with given identifier, default value, and lower and upper bound.
     * 
     * @param theid the identifier
     * @param thedefaultValue the default value
     * @param thelowerBound the lower bound, or {@code null} if the default lower bound shall be taken
     * @param theupperBound the upper bound
     */
    public Property(final String theid, final T thedefaultValue,
            final Comparable<? super T> thelowerBound, final Comparable<? super T> theupperBound) {
        this(theid, thedefaultValue, thelowerBound);
        if (theupperBound != null) {
            this.upperBound = theupperBound;
        }
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof IProperty<?>) {
            return this.id.equals(((IProperty<?>) obj).getId());
        } else {
            return false;
        }
    }
    
    @Override
    public int hashCode() {
        return id.hashCode();
    }
    
    @Override
    public String toString() {
        return id;
    }
    
    /**
     * Returns the default value of this property.
     * 
     * @return the default value.
     */
    @SuppressWarnings("unchecked")
    public T getDefault() {
        // Clone the default value if it's a Cloneable. One would have to use reflection for this to work
        // properly (classes implementing Cloneable are not required to make their clone() method
        // public, so we need to check if they have such a method and invoke it via reflection, which
        // results in ugly and unchecked type casting)
        // HOWEVER, since GWT doesn't support reflection, a helper class is used that is able to clone 
        // those objects that are viable as property values 
        if (defaultValue instanceof Cloneable) {
            Object clone = ElkReflect.clone(defaultValue);
            if (clone == null) {
                throw new IllegalStateException(
                        "Couldn't clone property '" + id + "'. " + "Make sure its type is registered with the "
                                + ElkReflect.class.getSimpleName() + " utility class.");
            }
            return (T) clone;
        } else {
            return defaultValue;
        }
    }
    
    @Override
    public String getId() {
        return id;
    }
    
    @Override
    public Comparable<? super T> getLowerBound() {
        return lowerBound;
    }

    @Override
    public Comparable<? super T> getUpperBound() {
        return upperBound;
    }

    @Override
    public int compareTo(final IProperty<?> other) {
        return id.compareTo((String) other.getId());
    }

    /** the default lower bound, which is smaller than everything else. */
    public static final Comparable<Object> NEGATIVE_INFINITY = new Comparable<Object>() {
        public int compareTo(final Object other) {
            // Ignore FindBugs warning
            return -1;
        }
        
        public String toString() {
            return "-\u221e";
        };
    };
    
    /** the default upper bound, which is greater than everything else. */
    public static final Comparable<Object> POSITIVE_INFINITY = new Comparable<Object>() {
        public int compareTo(final Object other) {
            // Ignore FindBugs warning
            return 1;
        }
        
        public String toString() {
            return "+\u221e";
        };
    };
}
