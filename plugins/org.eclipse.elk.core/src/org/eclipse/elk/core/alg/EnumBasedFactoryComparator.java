/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.alg;

import java.util.Comparator;

/**
 * Comparator which sorts {@link ILayoutProcessorFactory layout processor factories} based on the order of their
 * enumeration constants. This only works for factories which are indeed enumerations. Plus, the ordinals of the
 * enumeration constants must be ordered such that they reflect dependencies between the processors.
 *
 * @param <G>
 *            the type of graph the processor will process.
 */
public final class EnumBasedFactoryComparator<G> implements Comparator<ILayoutProcessorFactory<G>> {

    /* (non-Javadoc)
     * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
     */
    @Override
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public int compare(final ILayoutProcessorFactory<G> factory1, final ILayoutProcessorFactory<G> factory2) {
        // If the two factories are of the same type, and if that type is an enumeration, we know how
        // to compare them
        Class<?> factory1Class = factory1.getClass();
        Class<?> factory2Class = factory2.getClass();

        // elkjs-exclude-start
        if (!factory1Class.isEnum() || (factory1Class.getClass() != factory2Class.getClass()
                && factory1Class.getDeclaringClass() != factory2Class.getDeclaringClass())) {

            throw new UnsupportedOperationException("This comparator can only compare enumeration "
                    + "constants that are part of the same enumeration.");
        }
        // elkjs-exclude-end

        return ((Enum) factory1).compareTo((Enum) factory2);
    }

}
