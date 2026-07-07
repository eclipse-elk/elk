/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * This is a special data structure made to efficiently move nodes between layers and still be able to efficiently get
 * all nodes of a layer. Moreover, it is order preserving.
 * It consists of a {@link LinkedHashMap} with a list of nodes to have a ordered list for each layer id and
 * a [{@link HashMap} to get the layer of a node.
 */
public class BiLinkedHashMultiMap<K extends Comparable<K>, V> {
    
    private LinkedHashMap<K, LinkedList<V>> multiMapKeyToValues = new LinkedHashMap<>();
    private LinkedHashMap<V, K> hashMapValuesToKey = new LinkedHashMap<V, K>();
    
    /**
     * Adds all values to the respective key.
     * 
     * @param key The key.
     * @param values A set of values
     */
    public void putAll(final K key, final List<V> values) {
        for (V value : values) {
            put(key, value);
        }
    }
    
    /**
     * Adds a new value linked to the given key and update prior occurrences of the value.
     * 
     * @param key The key.
     * @param value The value.
     */
    public void put(final K key, final V value) {
        // Remove old value.
        K oldKey = hashMapValuesToKey.get(value);
        if (oldKey != null) {
            LinkedList<V> values = multiMapKeyToValues.get(oldKey);
            values.remove(value);
        }
        // Add new value
        LinkedList<V> values = multiMapKeyToValues.get(key);
        if (values == null) {
            values = new LinkedList<>();
            multiMapKeyToValues.put(key, values);
        }
        values.add(value);
        
        hashMapValuesToKey.put(value, key);
    }
    
    /**
     * Returns the key associated with the given value.
     * 
     * @param value The value.
     * @return The key.
     */
    public K getKey(final V value) {
        return hashMapValuesToKey.get(value);
    }
    
    /**
     * Returns the values associated with the given key.
     * 
     * @param key The key.
     * @return The set of values for the given key.
     */
    public LinkedList<V> getValues(final K key) {
        LinkedList<V> toReturn = multiMapKeyToValues.get(key);
        if (toReturn == null) {
            toReturn = new LinkedList<>();
        }
        return toReturn;
    }
    
    /**
     * Returns a set of all key.
     * 
     * @return The keys.
     */
    public Set<K> keySet() {
        return multiMapKeyToValues.keySet();
    }
    
    public boolean isMaximalKey(K key) {
        for (K otherKey : multiMapKeyToValues.keySet()) {
            if (key.compareTo(otherKey) < 0) {
                return false;
            }
        }
        return true;
    }
    
    public boolean isMinimalKey(K key) {
        for (K otherKey : multiMapKeyToValues.keySet()) {
            if (key.compareTo(otherKey) > 0) {
                return false;
            }
        }
        return true;
    }
    
    
}
