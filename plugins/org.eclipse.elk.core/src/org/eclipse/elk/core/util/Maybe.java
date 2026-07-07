/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Object that may contain another object, inspired by the Haskell type <i>Maybe</i>.
 * <p>
 * This class can be used to wrap objects for anonymous classes:
 * <pre>
 * Myclass foo() {
 *     final Maybe<Myclass> maybe = new Maybe<Myclass>();
 *     PlatformUI.getWorkbench().getDisplay().syncExec(new Runnable() {
 *         public void run() {
 *             maybe.set(new Myclass("bar"));
 *         }
 *     });
 *     return maybe.get();
 * }
 * </pre>
 * </p>
 * <p>
 * Another use is as a wrapper for synchronization on objects that may be {@code null}:
 * <pre>
 * Maybe<Myclass> maybe = new Maybe<Myclass>();
 * 
 * void thread1() {
 *     maybe.set(new Myclass("foo"));
 *     synchronized (maybe) {
 *         maybe.notify();
 *     )
 * }
 * 
 * void thread2() {
 *     synchronized (maybe) {
 *         if (maybe.get() == null) {
 *             maybe.wait();
 *         }
 *     }
 *     maybe.get().bar();
 * }
 * </pre>
 * </p>
 * 
 * @param <T> type of contained object
 * @author msp
 */
public final class Maybe<T> implements Iterable<T> {
    
    /**
     * Create a maybe with inferred generic type.
     * 
     * @param <D> the generic type
     * @return a new instance of given type
     */
    public static <D> Maybe<D> create() {
        return new Maybe<D>();
    }

    /** the contained object, which may be {@code null}. */
    private T object;
    
    /**
     * Creates a maybe without an object.
     */
    public Maybe() {
        this.object = null;
    }
    
    /**
     * Creates a maybe with the given object.
     * 
     * @param theobject the object to contain
     */
    public Maybe(final T theobject) {
        this.object = theobject;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof Maybe<?>) {
            Maybe<?> other = (Maybe<?>) obj;
            return this.object == null ? other.object == null
                    : this.object.equals(other.object);
        } else {
            return false;
        }
    }
    
    @Override
    public int hashCode() {
        if (object == null) {
            return 0;
        } else {
            return object.hashCode();
        }
    }
    
    @Override
    public String toString() {
        if (object == null) {
            return "maybe(null)";
        } else {
            return "maybe(" + object.toString() + ")";
        }
    }

    /**
     * Sets the contained object.
     *
     * @param theobject the object to set
     */
    public void set(final T theobject) {
        this.object = theobject;
    }

    /**
     * Returns the contained object.
     *
     * @return the contained object
     */
    public T get() {
        return object;
    }

    @Override
    public Iterator<T> iterator() {
        return new Iterator<T>() {
            private boolean visited = false;
            public boolean hasNext() {
                return !visited && object != null;
            }
            public T next() {
                if (visited || object == null) {
                    throw new NoSuchElementException();
                }
                visited = true;
                return object;
            }
            public void remove() {
                if (!visited || object == null) {
                    throw new IllegalStateException();
                }
                object = null;
            }
        };
    }

    /**
     * Clear any contained object.
     */
    public void clear() {
        object = null;
    }
    
    /**
     * Determine whether any object is contained.
     * 
     * @return false if an object instance is contained, true otherwise
     */
    public boolean isEmpty() {
        return object == null;
    }
    
}
