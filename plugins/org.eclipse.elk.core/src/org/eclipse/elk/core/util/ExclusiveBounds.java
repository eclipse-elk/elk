/*******************************************************************************
 * Copyright (c) 2016 TypeFox and others.
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
package org.eclipse.elk.core.util;

/**
 * Utility {@link Comparable} implementations that can be used to specify exclusive upper and lower
 * bounds for layout options. For instance, a layout option whose values may be in the range (0, 1)
 * could be specified with the lower bound {@code greaterThan(0)} and the upper bound {@code lessThan(1)}.
 * 
 * <p>For inclusive bounds, you can simply use the numbers themselves as lower and upper bounds.</p>
 */
public final class ExclusiveBounds {
    
    /**
     * Hidden default constructor to avoid instantiation.
     */
    private ExclusiveBounds() {
    }
    
    /**
     * Create a lower bound that does not include the limit.
     * 
     * @param exclusiveLowerBound the lower bound
     * @return a comparable that excludes the given lower bound
     */
    public static Comparable<Number> greaterThan(final double exclusiveLowerBound) {
        return new ExclusiveLowerBound(exclusiveLowerBound);
    }
    
    /**
     * Lower bound that does not include the limit.
     */
    public static class ExclusiveLowerBound implements Comparable<Number> {
        
        private final double exclusiveLowerBound;
        
        /**
         * Create an exclusive lower bound.
         */
        public ExclusiveLowerBound(final double exclusiveLowerBound) {
            this.exclusiveLowerBound = exclusiveLowerBound;
        }

        @Override
        public int compareTo(final Number x) {
            if (exclusiveLowerBound < x.doubleValue()) {
                return -1;
            } else {
                return 1;
            }
        }
        
        @Override
        public boolean equals(final Object obj) {
            if (obj instanceof ExclusiveLowerBound) {
                ExclusiveLowerBound other = (ExclusiveLowerBound) obj;
                return this.exclusiveLowerBound == other.exclusiveLowerBound;
            }
            return false;
        }
        
        @Override
        public int hashCode() {
            return Double.hashCode(exclusiveLowerBound);
        }
        
        @Override
        public String toString() {
            return exclusiveLowerBound + " (exclusive)";
        }
        
    }
    
    /**
     * Create an upper bound that does not include the limit.
     * 
     * @param exclusiveUpperBound the upper bound
     * @return a comparable that excludes the given upper bound
     */
    public static Comparable<Number> lessThan(final double exclusiveUpperBound) {
        return new ExclusiveUpperBound(exclusiveUpperBound);
    }
    
    /**
     * Upper bound that does not include the limit.
     */
    public static class ExclusiveUpperBound implements Comparable<Number> {
        
        private final double exclusiveUpperBound;
        
        /**
         * Create an exclusive upper bound.
         */
        public ExclusiveUpperBound(final double exclusiveUpperBound) {
            this.exclusiveUpperBound = exclusiveUpperBound;
        }

        @Override
        public int compareTo(final Number x) {
            if (exclusiveUpperBound > x.doubleValue()) {
                return 1;
            } else {
                return -1;
            }
        }
        
        @Override
        public boolean equals(final Object obj) {
            if (obj instanceof ExclusiveUpperBound) {
                ExclusiveUpperBound other = (ExclusiveUpperBound) obj;
                return this.exclusiveUpperBound == other.exclusiveUpperBound;
            }
            return false;
        }
        
        @Override
        public int hashCode() {
            return Double.hashCode(exclusiveUpperBound);
        }
        
        @Override
        public String toString() {
            return exclusiveUpperBound + " (exclusive)";
        }
        
    }

}
