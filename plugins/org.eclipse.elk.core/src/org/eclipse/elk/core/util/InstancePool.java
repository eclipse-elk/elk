/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.LinkedList;

/**
 * A pool for class instances. The pool can hold a configurable number of instances of the class. Note
 * that the number of instances created can well exceed the capacity of this pool. The capacity only
 * influences how many of those instances can be kept for later reuse.
 * 
 * <p>All of the instance methods in this class are thread-safe.</p>
 *
 * @param <T> the type of instances that are held by this pool
 * @author msp
 */
public final class InstancePool<T> {
    
    /** an infinite amount of instances. */
    public static final int INFINITE = -1;
    
    /** the instance factory to use for this pool. */
    private final IFactory<T> factory;
    /** the list of currently held instances. */
    private final LinkedList<T> instances = new LinkedList<T>();
    /** the configured instance limit. */
    private int limit;
    /** lock used for synchronization. */
    private final Object lock = new Object();
    
    /**
     * Create an instance pool with an infinite capacity. 
     * 
     * @param thefactory the instance factory
     */
    public InstancePool(final IFactory<T> thefactory) {
        this(thefactory, INFINITE);
    }
    
    /**
     * Create an instance pool with given capacity.
     * 
     * @param thefactory the instance factory
     * @param thelimit the maximal number of instances that shall be kept in the pool
     */
    public InstancePool(final IFactory<T> thefactory, final int thelimit) {
        assert thefactory != null;
        this.factory = thefactory;
        this.limit = thelimit;
    }
    
    /**
     * Fetch an instance from the pool. If no instance is available, a new one is created.
     * 
     * @return a class instance
     */
    public T fetch() {
        synchronized (lock) {
            if (instances.isEmpty()) {
                return factory.create();
            }
            return instances.removeFirst();
        }
    }
    
    /**
     * Release an instance into the pool to be used again unless the pool's capacity is already reached.
     * Only instances that are still usable may be released.
     * 
     * @param obj a class instance
     */
    public void release(final T obj) {
        synchronized (lock) {
            if (limit < 0 || instances.size() < limit) {
                instances.addLast(obj);
            } else {
                factory.destroy(obj);
            }
        }
    }
    
    /**
     * Clear the instance pool by disposing all instances that are currently held.
     */
    public void clear() {
        synchronized (lock) {
            for (T obj : instances) {
                factory.destroy(obj);
            }
            instances.clear();
        }
    }

}
