/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.util;

import java.lang.annotation.Annotation;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.List;

import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * This class contains utility methods useful for other classes of the test framework.
 */
public final class TestUtil {

    /**
     * Private constructor to avoid instantiation.
     */
    private TestUtil() {
    }

    /**
     * Adds an exception to the list of errors if the given method is not public.
     */
    public static boolean ensurePublic(final FrameworkMethod method, final List<Throwable> errors) {
        if (!method.isPublic()) {
            errors.add(new Exception("Method " + method.getName() + " must be public"));
            return false;
        } else {
            return true;
        }
    }

    /**
     * Adds an exception to the list of errors if the given method has a non-void return type.
     */
    public static boolean ensureReturnsType(final FrameworkMethod method, final List<Throwable> errors,
            final Class<?> expectedType) {

        if (!expectedType.isAssignableFrom(method.getReturnType())) {
            errors.add(new Exception(
                    "Method " + method.getName() + " must have return type " + expectedType.getCanonicalName()));
            return false;
        } else {
            return true;
        }
    }

    /**
     * Adds an exception to the list of errors if the given method does not have parameters that can accept arguments of
     * the given types.
     */
    public static boolean ensureParameters(final FrameworkMethod method, final List<Throwable> errors,
            final Class<?>... expectedTypes) {

        Method javaMethod = method.getMethod();
        boolean checkPassed = true;

        if (javaMethod.getParameterCount() == expectedTypes.length) {
            Class<?>[] parameterTypes = javaMethod.getParameterTypes();

            for (int i = 0; i < expectedTypes.length; i++) {
                if (!parameterTypes[i].isAssignableFrom(expectedTypes[i])) {
                    checkPassed = false;
                }
            }

        } else {
            checkPassed = false;
        }

        if (checkPassed) {
            return true;
        } else {
            errors.add(new Exception("Method " + method.getName() + " must expect parameters of types: "
                    + Arrays.toString(expectedTypes)));
            return false;
        }
    }

    /**
     * Adds exceptions to the list of errors if the given method has annotations of the given annotation types.
     */
    @SafeVarargs
    public static boolean ensureNotAnnotatedWith(final FrameworkMethod method, final List<Throwable> errors,
            final Class<? extends Annotation>... forbidden) {

        boolean success = true;

        for (Annotation annotation : method.getAnnotations()) {
            for (Class<? extends Annotation> forbiddenAnnotationType : forbidden) {
                if (forbiddenAnnotationType.equals(annotation.getClass())) {
                    errors.add(new Exception("Method " + method.getName() + " cannot have annotation "
                            + forbiddenAnnotationType.getName()));
                    success = false;
                }
            }
        }

        return success;
    }

    /**
     * Adds an exception to the list of errors if the given class does not have a single public no-parameters
     * constructor.
     */
    public static boolean ensureSinglePublicNoParameterConstructor(final TestClass testClass,
            final List<Throwable> errors) {

        Class<?> javaClass = testClass.getJavaClass();
        Constructor<?>[] constructors = javaClass.getConstructors();

        if (constructors.length != 1 || constructors[0].getParameterCount() != 0) {
            errors.add(new Exception("Test class must have a single public constructor without parameters."));
            return false;
        }

        return true;
    }

    /**
     * Tries to create and return an instance of the test class.
     * 
     * @return an instance of the test class.
     * @throws Throwable
     *             if anything goes wrong.
     */
    public static Object createTestClassInstance(final TestClass testClass) throws Throwable {
        return testClass.getOnlyConstructor().newInstance();
    }

}
