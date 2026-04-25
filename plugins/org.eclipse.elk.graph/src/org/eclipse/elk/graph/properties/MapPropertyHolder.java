/*******************************************************************************
 * Copyright (c) 2009, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.properties;

import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * An implementation of {@link IPropertyHolder} based on a {@link HashMap}.
 *
 * @author msp
 */
public class MapPropertyHolder implements IPropertyHolder, Serializable {

    /** the serial version UID. */
    private static final long serialVersionUID = 4507851447415709893L;
    
    /** map of property identifiers to their values. */
    private HashMap<IProperty<?>, Object> propertyMap;
    
    @Override
    public <T> MapPropertyHolder setProperty(final IProperty<? super T> property, final T value) {
        if (value == null) {
            getProperties().remove(property);
        } else {
            getProperties().put(property, value);
        }
        
        return this;
    }
    
    @Override
    @SuppressWarnings("unchecked")
    public <T> T getProperty(final IProperty<T> property) {
        Object value = getProperties().get(property);
        if (value instanceof IPropertyValueProxy) {
            value = ((IPropertyValueProxy) value).resolveValue(property);
            if (value != null) {
                getProperties().put(property, value);
                return (T) value;
            }
        } else if (value != null) {
            return (T) value;
        }

        // the reason for the side effect below is that if a default value has been returned 
        // and the object is altered by the user, the user expects the altered object to be 
        // the value of the property in case he asks for the property again
        
        // Retrieve the default value and memorize it for our property
        T defaultValue = property.getDefault();
        if (defaultValue instanceof Cloneable) {
            // We are now dealing with a clone of the default value which me may safely store away
            // for further modification
            setProperty(property, defaultValue);
        }
        return defaultValue;
    }
    
    @Override
    public boolean hasProperty(final IProperty<?> property) {
        return propertyMap != null && propertyMap.containsKey(property);
    }
    
    @Override
    public MapPropertyHolder copyProperties(final IPropertyHolder other) {
        if (other == null) {
            return this;
        }

        final Map<IProperty<?>, Object> otherMap = other.getAllProperties();
        if (!otherMap.isEmpty()) {
            if (this.propertyMap == null) {
                propertyMap = new HashMap<IProperty<?>, Object>(otherMap);
            } else {
                this.propertyMap.putAll(otherMap);
            }
        }

        return this;
    }
    
    @Override
    public Map<IProperty<?>, Object> getAllProperties() {
        if (propertyMap == null) {
            return Collections.emptyMap();
        } else {
            return propertyMap;
        }
    }
    
    /**
     * Returns the property map, creating a new map if there hasn't been one so far.
     * 
     * @return the property map.
     */
    private Map<IProperty<?>, Object> getProperties() {
        if (propertyMap == null) {
            propertyMap = new HashMap<>();
        }
        return propertyMap;
    }
    
}
