/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.util;

import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Provides the means to register a new instance function and a clone function with class types. The functionality to
 * create new instances and clone objects of certain types is required by ELK's property mechanism. The rather
 * cumbersome solution to use this class is mandatory because we want to stay compatible with GWT which does not support
 * Java's reflection and clonable mechanisms. Otherwise using {@link Class#newInstance()} and retrieving a class's
 * {@code clone()} method via reflection would be enough.
 * 
 * <p>
 * Special treatment is devoted to {@link EnumSet}s. To avoid having to register the concrete implementations of
 * {@link EnumSet}.
 * </p>
 * 
 * <h3>Notes</h3>
 * <p>
 * The class has to reside in the elk.graph plugin since it is used by the
 * {@link org.eclipse.elk.graph.properties.Property Property} class. The majority of necessary classes are 
 * registered in the constructor of ELK's {@code LayoutMetaDataService}. 
 * </p>
 */
public final class ElkReflect {

    private static final Map<Class<?>, NewInstanceFunction> REGISTRY_NEW = new HashMap<>();
    private static final Map<Class<?>, CloneFunction> REGISTRY_CLONE = new HashMap<>();
    
    static {
        // special treatment and no new instance function required for EnumSets, register clone here
        ElkReflect.registerClone(EnumSet.class, (es) -> ((EnumSet<?>) es).clone());
    }
    
    private ElkReflect() { }

    /**
     * Register both a new instance function and a clone function for the passed {@code clazz}.
     */
    public static void register(final Class<?> clazz, final NewInstanceFunction newFun, final CloneFunction cloneFun) {
        if (clazz != null && newFun != null) {
            REGISTRY_NEW.put(clazz, newFun);
        }
        if (clazz != null && cloneFun != null) {
            REGISTRY_CLONE.put(clazz, cloneFun);
        }
    }
    
    /**
     * Register a new instance function for the passed {@code clazz}.
     */
    public static void registerNewInstance(final Class<?> clazz, final NewInstanceFunction newFun) {
        register(clazz, newFun, null);
    }
    
    /**
     * Register a clone function for the passed {@code clazz}.
     * @param cloneFun
     */
    public static void registerClone(final Class<?> clazz, final CloneFunction cloneFun) {
        register(clazz, null, cloneFun);
    }
    
    
    /**
     * If a {@link NewInstanceFunction} is registered for {@code clazz}, the function is called and the created instance
     * is returned.
     * 
     * @return either a new instance of the requested class type or {@code null} if the class is unknown.
     */
    @SuppressWarnings("unchecked")
    public static <T> T newInstance(final Class<T> clazz) {
        if (REGISTRY_NEW.containsKey(clazz)) {
            return (T) REGISTRY_NEW.get(clazz).newInstance();
        } else {
            return null;
        }
    }
    
    /**
     * Clones the passed {@code clonee} if a {@link CloneFunction} is registered for {@code clonee}'s class. 
     * Note that for {@link EnumSet}s an {@code instanceof} check is performed, for all other classes the 
     * exact type has to be known. 
     * 
     * @param clonee
     *            the object to clone
     * @return either a cloned version of {@code clonee} or {@code null} if {@code clonee}'s class is unknown.
     */
    @SuppressWarnings("unchecked")
    public static <T> T clone(final T clonee) {
        // special treatment for EnumSets
        if (clonee instanceof EnumSet<?>) {
            return (T) REGISTRY_CLONE.get(EnumSet.class).clone(clonee);
        }
        
        if (REGISTRY_CLONE.containsKey(clonee.getClass())) {
            return (T) REGISTRY_CLONE.get(clonee.getClass()).clone(clonee);
        }
        return null;
    }
    
    /**
     * Functional interface to specify a function that knows how to create a new instance of a particular class.
     */
    @FunctionalInterface
    public interface NewInstanceFunction {
        /** Create and return a new instance of a specific class. */
        Object newInstance();
    }
    
    /**
     * Functional interface to specify a function that knows how to clone a particular object.  
     */
    @FunctionalInterface
    public interface CloneFunction {
        /** Create a cloned version of the passed object. */
        Object clone(Object o);
    }
}
