/*******************************************************************************
 * Copyright (c) 2009, 2017 Kiel University and others.
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
package org.eclipse.elk.graph.properties;

import java.util.Map;

/**
 * Interface for holders of property values.
 *
 * @author msp
 */
public interface IPropertyHolder {
    
    /**
     * Sets a property value. No type checking is performed while setting, so
     * users of this method must take care that the right object types are generated.
     * 
     * @param <T> type of property
     * @param property the property to set
     * @param value the new value
     * @return <code>this</code> {@link IPropertyHolder} for convenience
     */
    <T> IPropertyHolder setProperty(IProperty<? super T> property, T value);
    
    /**
     * Retrieves a property value. If the property is not set, its default value shall be taken,
     * which is taken from the given property instance.
     * 
     * @param <T> type of property
     * @param property the property to get
     * @return the current value, or the default value if the property is not set
     */
    <T> T getProperty(IProperty<T> property);
    
    /**
     * Checks whether a value is configured for the given property. If not, the next call to
     * {@link #getProperty(IProperty)} will return the property's default value and set the
     * property's value to the default value for this property holder, provided that the
     * default value is not {@code null}. After that, the property value is set and thus
     * all further calls to {@link #hasProperty(IProperty)} will return {@code true} for
     * that property.
     * 
     * @param property the property.
     * @return {@code true} or {@code false} as a value is or is not set for the property. 
     */
    boolean hasProperty(IProperty<?> property);
    
    /**
     * Copy all properties from another property holder to this one.
     * 
     * @param holder another property holder
     * @return <code>this</code> {@link IPropertyHolder} for convenience
     */
    IPropertyHolder copyProperties(IPropertyHolder holder);
    
    /**
     * Returns a map of all assigned properties with associated values.
     * 
     * @return a map of all properties
     */
    Map<IProperty<?>, Object> getAllProperties();

}
